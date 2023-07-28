import sys
sys.path.append("/home/david/catkin_ws/src/franka_david/scripts/qdc-manip")

import time
import numpy as np
import cv2
from scipy import optimize

from kinect import KinectClient
from franka import Franka, DEFAULT_ORN


# Default ORN is DEFAULT_ORN = [2.22, 2.22, 0.0]
measured_pts, observed_pts, observed_pix, world2camera = [None]*4


def calibrate(cam, franka, workspace_bounds, ee_to_checker=0.36, calib_grid_step=0.1):
    global measured_pts, observed_pts, observed_pix, world2camera
    # Constants
    # Distance to checker includes the distance to the fingers
    ee_to_checker = 0.355
    checkerboard_offset = np.array([0, 0, ee_to_checker])
    sleep_time = 5

    # Construct 3D calibration grid across workspace
    gridspace_x = np.linspace(workspace_bounds[0, 0],
                              workspace_bounds[0, 1],
                              1 + int((workspace_bounds[0, 1]-workspace_bounds[0, 0])/calib_grid_step))
    gridspace_y = np.linspace(workspace_bounds[1, 0],
                              workspace_bounds[1, 1],
                              1+int((workspace_bounds[1, 1]-workspace_bounds[1, 0])/calib_grid_step))

    calib_grid_x, calib_grid_y, calib_grid_z = np.meshgrid(gridspace_x, gridspace_y, workspace_bounds[2, 0])

    num_calib_grid_pts = calib_grid_x.shape[0] * calib_grid_x.shape[1]*calib_grid_x.shape[2]
    calib_grid_x.shape = (num_calib_grid_pts, 1)
    calib_grid_y.shape = (num_calib_grid_pts, 1)
    calib_grid_z.shape = (num_calib_grid_pts, 1)
    calib_grid_pts = np.concatenate(
        (calib_grid_x, calib_grid_y, calib_grid_z), axis=1)

    # Move robot to each calibration point in workspace
    measured_pts = list()
    observed_pts = list()
    observed_pix = list()
    prev_point_y = 0.0
    print(f"The number of calibgration grid points is {num_calib_grid_pts}")
    for calib_pt_idx in range(num_calib_grid_pts):
        tool_position = calib_grid_pts[calib_pt_idx, :]
        tool_position[2] = workspace_bounds[2, 1]
        print(f"Moving to position={tool_position}")
        new_point_y = tool_position[1]

        if prev_point_y != new_point_y:
            prev_point_y = new_point_y
            duration = 5
            franka.send_home()
            time.sleep(2)
        else:
            duration = sleep_time
        franka.movel(params=[list(tool_position) + DEFAULT_ORN],
                     traj_duration=duration)

        time.sleep(duration)

        while True:
            print("Trying to find checker")
            color_im, depth_im = cam.get_rgbd(repeats=10)
            chckr_size = (3, 3)
            refine_criteria = (cv2.TERM_CRITERIA_EPS +
                               cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
            bgr_im = cv2.cvtColor(color_im, cv2.COLOR_RGB2BGR)
            gray_im = cv2.cvtColor(bgr_im, cv2.COLOR_RGB2GRAY)
            chckr_found, crnrs = cv2.findChessboardCorners(
                gray_im, chckr_size,
                # None, 0
                None, cv2.CALIB_CB_ADAPTIVE_THRESH
            )
            if chckr_found:
                crnrs_refined = cv2.cornerSubPix(
                    gray_im, crnrs, (3, 3), (-1, -1), refine_criteria)
                block_pix = crnrs_refined[4, 0, :]
                print("Checker found")
                break
            time.sleep(0.01)

        # Get observed checkerboard center 3D point in camera space
        block_z = depth_im[int(np.round(block_pix[1])), int(np.round(block_pix[0]))]
        block_x = np.multiply(
            block_pix[1] - cam.color_intr[0, 2],
            block_z / cam.color_intr[0, 0]
        )
        block_y = np.multiply(
            block_pix[0] - cam.color_intr[1, 2],
            block_z / cam.color_intr[1, 1]
        )
        if block_z == 0:
            continue

        # Save calibration point and observed checkerboard center
        observed_pts.append([block_x, block_y, block_z])
        tool_position += checkerboard_offset
        measured_pts.append(tool_position)
        observed_pix.append(block_pix)

        # Draw and display the corners
        center = np.array(block_pix).astype(np.int16)
        # print(f"Center is={center}, type={type(center)}")
        vis_im = cv2.circle(color_im, tuple(center), 7, (0, 255, 0), 2)
        cv2.imshow('Calibration', cv2.cvtColor(vis_im, cv2.COLOR_RGB2BGR))
        cv2.waitKey(2)

    # Move robot back to home pose
    # ur5.homej(blocking=True)

    measured_pts = np.asarray(measured_pts)
    observed_pts = np.asarray(observed_pts)
    observed_pix = np.asarray(observed_pix)
    # Save the arrays in case we get an error
    np.save('measured_pts.npy', measured_pts)
    np.save('observed_pts.npy', observed_pts)
    np.save('observed_pix.npy', observed_pix)
    world2camera = np.eye(4)

    # Estimate rigid transform with SVD (from Nghia Ho)
    def get_rigid_transform(A, B):
        assert len(A) == len(B)
        N = A.shape[0]  # Total points
        centroid_A = np.mean(A, axis=0)
        centroid_B = np.mean(B, axis=0)
        AA = A - np.tile(centroid_A, (N, 1))  # Centre the points
        BB = B - np.tile(centroid_B, (N, 1))
        # Dot is matrix multiplication for array
        H = np.dot(np.transpose(AA), BB)
        U, S, Vt = np.linalg.svd(H)
        R = np.dot(Vt.T, U.T)
        if np.linalg.det(R) < 0:  # Special reflection case
            Vt[2, :] *= -1
            R = np.dot(Vt.T, U.T)
        t = np.dot(-R, centroid_A.T) + centroid_B.T
        return R, t

    def get_rigid_transform_error(z_scale):
        global measured_pts, observed_pts, observed_pix, world2camera

        # Apply z offset and compute new observed points
        # using camera intrinsics
        observed_z = observed_pts[:, 2:] * z_scale
        observed_x = np.multiply(
            observed_pix[:, [0]]-cam.color_intr[0, 2],
            observed_z/cam.color_intr[0, 0])
        observed_y = np.multiply(
            observed_pix[:, [1]]-cam.color_intr[1, 2],
            observed_z/cam.color_intr[1, 1])
        new_observed_pts = np.concatenate(
            (observed_x, observed_y, observed_z), axis=1)

        # Estimate rigid transform between measured points
        # and new observed points
        R, t = get_rigid_transform(np.asarray(
            measured_pts), np.asarray(new_observed_pts))
        t.shape = (3, 1)
        world2camera = np.concatenate(
            (np.concatenate((R, t), axis=1), np.array([[0, 0, 0, 1]])), axis=0)

        # Compute rigid transform error
        registered_pts = np.dot(R, np.transpose(
            measured_pts)) + np.tile(t, (1, measured_pts.shape[0]))
        error = np.transpose(registered_pts) - new_observed_pts
        error = np.sum(np.multiply(error, error))
        rmse = np.sqrt(error/measured_pts.shape[0])
        return rmse

    # Optimize z scale w.r.t. rigid transform error
    print('Calibrating...')
    z_scale_init = 1
    optim_result = optimize.minimize(
        get_rigid_transform_error,
        np.asarray(z_scale_init),
        method='Nelder-Mead')
    camera_depth_offset = optim_result.x

    # Save camera optimized offset and camera pose
    print('Saving calibration files...')
    np.savetxt('camera_depth_scale.txt',
               camera_depth_offset,
               delimiter=' ')
    get_rigid_transform_error(camera_depth_offset)
    camera_pose = np.linalg.inv(world2camera)
    return camera_pose


if __name__ == "__main__":

    workspace_bounds = np.array([
        [0.30053, 0.72368],
        # [-0.35026, 0.35026],
        [-0.35026, 0.10026],
        [0.13, 0.13]
    ])
    # ur5_left = get_single_ur5(ur5_ip="192.168.56.101")
    franka = Franka(init_node=True)
    franka.send_home()


    max_depth = 2.0  # TODO To change to actual maximum depth to be used
    gundam_ip = "130.233.123.198"

    cam = KinectClient(rotate_img=False, ip=gundam_ip)

    print(f"Camera intrinsics={cam.color_intr}")

    np.savetxt('top_down_left_ur5_cam_pose.txt',
               calibrate(cam, franka, workspace_bounds),
               delimiter=' ')

