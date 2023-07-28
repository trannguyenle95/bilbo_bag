# Copyright 2021 (c) Aalto University - All Rights Reserved
# Author: David Blanco Mulero <david.blancomulero@aalto.fi>
#

import numpy as np
from scipy.spatial.transform import Rotation as R


def check_single_action_reachability(arm_base: np.array, p1: np.array, reach_limit: np.array):
    return check_arm_reachability(reach_limit, arm_base, p1)


def check_action_reachability(action: str, arm_base: np.array,
                              p1: np.array, p2: np.array, reach_limit: np.array):
    if action == 'fling':
        # right and left must reach each point respectively
        return check_arm_reachability(reach_limit, arm_base, p1), None
    elif action == 'drag' or action == 'place' or action == 'dynamic' or action == 'fling-place':
        # either right can reach both or left can reach both
        if check_arm_reachability(reach_limit, arm_base, p1) and \
                check_arm_reachability(reach_limit, arm_base, p2):
            return True, 'left'
        else:
            return False, None
    raise NotImplementedError()


def check_arm_reachability(reach_limit: np.array, arm_base: np.array, reach_pos: np.array):
    return np.linalg.norm(arm_base - reach_pos) < reach_limit


def get_action_params(action_primitive: str, x: int, y: int, z: int,
                      pix_grasp_dist: int, pix_drag_dist: int, pix_place_dist: int):
    p1 = np.array([y, z])
    p2 = p1.copy()
    if action_primitive == 'fling':
        p1[0] = p1[0] + pix_grasp_dist
        p2[0] = p2[0] - pix_grasp_dist
    elif action_primitive == 'drag':
        p2[0] += pix_drag_dist
    elif action_primitive == 'place' or action_primitive == 'fling-place' or action_primitive == 'dynamic':
        p2[0] += pix_place_dist
    else:
        raise Exception(f'Action Primitive not supported: {action_primitive}')
    return np.array([p1, p2])


def get_quintic_poly_coeffs(q0, qf, tf):
    v0, vf, a0, af = 0, 0, 0, 0
    known_vec = np.array([q0, v0, a0, qf, vf, af])
    M = np.array([[1, 0, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0, 0],
                  [0, 0, 2, 0, 0, 0],
                  [1, tf, tf ** 2, tf ** 3, tf ** 4, tf ** 5],
                  [0, 1, 2 * tf, 3 * tf ** 2, 4 * tf ** 3, 5 * tf ** 4],
                  [0, 0, 2, 6 * tf, 12 * tf ** 2, 20 * tf ** 3]])
    return np.dot(np.linalg.inv(M), known_vec)


def get_xyz_coordinates_given_tf(x0, xf, z0, zf, tf, dt):
    dx = (xf - x0) / 2
    dz = (zf - z0) / 2
    r = np.linalg.norm(np.array([dx, dz]))

    coefs = get_quintic_poly_coeffs(q0=np.pi, qf=0, tf=tf)
    t = np.arange(start=0, step=dt, stop=tf + dt)

    theta = coefs[0] + coefs[1] * t + coefs[2] * t ** 2 + coefs[3] * t ** 3 + coefs[4] * t ** 4 + coefs[5] * t ** 5
    xnew = x0 + dx + (dx * np.cos(theta))  # dx/dt =  dx/dtheta * dtheta/dt || np.sin(theta) * (coeff[1] * t  + ...)
    znew = z0 + dz + (dz * np.cos(theta))
    # xnew = x0 + (dx * np.cos(theta))
    # znew = z0 + (dz * np.cos(theta))
    ynew = r * np.sin(theta)
    # print("coefs: \n", coefs)
    return xnew, ynew, znew, theta, coefs


def get_mid_theta_vel(coefs, tf):
    t = tf / 2
    return coefs[1] + 2 * coefs[2] * t + 3 * coefs[3] * t ** 2 + 4 * coefs[4] * t ** 3 + 5 * coefs[5] * t ** 4


def get_xyz_coordinates_given_mid_theta_ref(start_pos, end_pos, mid_theta_ref):
    x0 = start_pos[0, 0]
    xf = end_pos[0, 0]
    z0 = start_pos[0, 2]
    zf = end_pos[0, 2]
    dt = 1 / 240

    tf = 0.01
    K = 0.005
    tol = 0.01
    it = 0
    dx = (xf - x0) / 2
    dz = (zf - z0) / 2
    r = np.linalg.norm(np.array([dx, dz]))
    if r <= 0.05:
        mid_theta_ref = -0.3

    x, y, z, theta, coefs = get_xyz_coordinates_given_tf(x0, xf, z0, zf, tf, dt)
    mid_theta = r * get_mid_theta_vel(coefs, tf)
    err = mid_theta_ref - mid_theta

    while (np.abs(err) > tol):
        tf += K * err
        x, y, z, theta, coefs = get_xyz_coordinates_given_tf(x0, xf, z0, zf, tf, dt)
        mid_theta = r * get_mid_theta_vel(coefs, tf)
        err = mid_theta_ref - mid_theta
        # print("iter: ", it, "mid_theta_ref: ", mid_theta_ref, "  mid_theta: ", mid_theta)
        it += 1
        if it > 10000:
            print(f"Used tf: {tf} sec. at iteration {it} for distance r: {r}.")
            break
    quintic_traj = np.array([x, y, z]).T

    return quintic_traj, theta, coefs, tf


def get_xyz_coordinates_given_mid_theta_ref_real(start_pos, end_pos, mid_theta_ref):
    x0 = start_pos[0, 0]
    xf = end_pos[0, 0]
    z0 = start_pos[0, 2]
    zf = end_pos[0, 2]
    dt = 1 / 1000  # Real should be 1000Hz  -- tf=1s, length trajectory=1200,

    tf = 10.0  # Initial estimate
    K = 0.05  # 0.005
    tol = 0.01
    it = 0
    dx = (xf - x0) / 2
    dz = (zf - z0) / 2
    r = np.linalg.norm(np.array([dx, dz]))
    if mid_theta_ref < -0.5:
        print("CAREFUL! THIS WILL POTENTIALLY SEND A REALLY FAST MOTION")
    if mid_theta_ref < -1.5:
        raise ValueError("Motion request is too fast")
    # if r <= 0.05:  # TODO Modify this value
    #     mid_theta_ref = -0.3

    x, y, z, theta, coefs = get_xyz_coordinates_given_tf(x0, xf, z0, zf, tf, dt)
    mid_theta = r * get_mid_theta_vel(coefs, tf)
    err = mid_theta_ref - mid_theta

    while (np.abs(err) > tol):
        tf += K * err
        x, y, z, theta, coefs = get_xyz_coordinates_given_tf(x0, xf, z0, zf, tf, dt)
        mid_theta = r * get_mid_theta_vel(coefs, tf)
        err = mid_theta_ref - mid_theta
        # print("iter: ", it, "mid_theta_ref: ", mid_theta_ref, "  mid_theta: ", mid_theta)
        it += 1
        if it > 10000:
            print(f"Used tf: {tf} sec. at iteration {it} for distance r: {r}.")
            break
    print(f"Found tf: {tf} sec. at iteration {it} for distance r: {r}.")
    quintic_traj = np.array([x, y, z]).T

    return quintic_traj, theta, coefs, tf


def get_xz_w_initial_conds(x0, xf, z0, zf, t0, tf, v0, vf, a0, af, dt):

    x_vec = np.array([x0, v0, a0, xf, vf, af])
    z_vec = np.array([z0, v0, a0, zf, vf, af])
    M = np.array([[1, t0, t0 ** 2, t0 ** 3, t0 ** 4, t0 ** 5],
                  [0, 1, 2 * t0, 3 * t0 ** 2, 4 * t0 ** 3, 5 * t0 ** 4],
                  [0, 0, 2, 6 * t0, 12 * t0 ** 2, 20 * t0 ** 3],
                  [1, tf, tf ** 2, tf ** 3, tf ** 4, tf ** 5],
                  [0, 1, 2 * tf, 3 * tf ** 2, 4 * tf ** 3, 5 * tf ** 4],
                  [0, 0, 2, 6 * tf, 12 * tf ** 2, 20 * tf ** 3]])

    x_coeffs, _, _, _ = np.linalg.lstsq(M, x_vec, rcond=None)
    z_coeffs, _, _, _ = np.linalg.lstsq(M, z_vec, rcond=None)

    t = np.arange(start=t0, step=dt, stop=tf + dt)
    xnew = x_coeffs[0] + x_coeffs[1] * (t - t0) + x_coeffs[2] * (t - t0) ** 2 + x_coeffs[3] * (t - t0) ** 3 + \
           x_coeffs[4] * (t - t0) ** 4 + x_coeffs[5] * (t - t0) ** 5
    znew = z_coeffs[0] + z_coeffs[1] * (t - t0) + z_coeffs[2] * (t - t0) ** 2 + z_coeffs[3] * (t - t0) ** 3 + \
           z_coeffs[4] * (t - t0) ** 4 + z_coeffs[5] * (t - t0) ** 5

    return xnew, znew, x_coeffs, z_coeffs


def compute_x(pos_start, pos_back, pos_front, pos_center, pos_last, tf1, tf2, tf3, tf4, dt):
    x0 = pos_start[0][0]
    xf = pos_back[0][0]
    z0 = pos_start[0][2]
    zf = pos_back[0][2]

    if pos_back > pos_start:
        sign_change = 1
    else:
        sign_change = -1

    t0 = 0
    v0 = 0
    a0 = 0
    vf = sign_change*0.7
    af = -sign_change*3.0

    xnew, _, _, _ = get_xz_w_initial_conds(x0, xf, z0, zf, t0, tf1, v0, vf, a0, af, dt)

    x0_2 = xf
    xf_2 = pos_front[0][0]
    z0_2 = zf
    zf_2 = pos_front[0][2]
    t0 = 0.0
    v0_2 = vf
    a0_2 = af
    vf_2 = sign_change*0.0
    af2 = sign_change*1.0

    xnew2, znew2, _, _ = get_xz_w_initial_conds(x0_2, xf_2, z0_2, zf_2, t0, tf2, v0_2, vf_2, a0_2, af2, dt)
    znew2[:] = zf_2

    x0_3 = xf_2
    xf_3 = pos_center[0][0]
    z0_3 = zf_2
    zf_3 = pos_center[0][2]
    t0 = 0.0
    v0_3 = vf_2
    a0_3 = af2
    vf_3 = sign_change*0.3
    # af3 = sign_change*1.0
    af3 = sign_change*0.5

    xnew3, _, _, _ = get_xz_w_initial_conds(x0_3, xf_3, z0_3, zf_3, t0, tf3, v0_3, vf_3, a0_3, af3, dt)

    x0_4 = xf_3
    xf_4 = pos_last[0][0]
    z0_4 = zf_3
    zf_4 = pos_last[0][2]
    t0 = 0.0
    v0_4 = vf_3
    a0_4 = af3
    vf_4 = 0.0
    af4 = 0.0

    xnew4, _, _, _ = get_xz_w_initial_conds(x0_4, xf_4, z0_4, zf_4, t0, tf4, v0_4, vf_4, a0_4, af4, dt)

    return xnew, xnew2, xnew3, xnew4


def compute_z(pos_start, pos_back, pos_front, pos_center, pos_last, tf1, tf3, tf4, dt):
    z0 = pos_start[0][2]
    zf = pos_back[0][2]
    zf_2 = pos_front[0][2]
    z0_3 = zf_2
    zf_3 = pos_center[0][2]

    t0 = 0
    v0 = 0
    a0 = 0
    vf = 0.0
    af = 0.0

    znew, _, _, _ = get_xz_w_initial_conds(z0, zf, z0, zf, t0, tf1, v0, vf, a0, af, dt)

    vf_2 = 0.0
    af2 = 0.0
    v0_3 = vf_2
    a0_3 = af2
    vf_3 = 0.0
    af3 = 0.02

    znew3, _, _, _ = get_xz_w_initial_conds(z0_3, zf_3, z0_3, zf_3, t0, tf3, v0_3, vf_3, a0_3, af3, dt)

    z0_4 = zf_3
    zf_4 = pos_last[0][2]
    v0_4 = vf_3
    a0_4 = af3
    vf_4 = 0.0
    af4 = 0.0

    znew4, _, _, _ = get_xz_w_initial_conds(z0_4, zf_4, z0_4, zf_4, t0, tf4, v0_4, vf_4, a0_4, af4, dt)
    znew4[:] = zf_4

    return znew, znew3, znew4


def compute_angle(angle, tf1, tf2, tf3, dt):
    p0 = angle[0]
    p1 = angle[1]
    a0 = 0
    t0 = 0.0
    v0 = 0.0
    vf = 10.8
    af = 180.0

    pitch_1, _, _, _ = get_xz_w_initial_conds(p0, p1, p0, p0, t0, tf1, v0, vf, a0, af, dt)

    p2 = angle[2]
    t0 = 0.0
    v0_2 = vf
    a0_2 = af
    vf_2 = -10.0
    af2 = -180.0

    pitch_2, _, _, _ = get_xz_w_initial_conds(p1, p2, p0, p0, t0, tf2, v0_2, vf_2, a0_2, af2, dt)

    p3 = angle[3]
    t0 = 0.0
    v0_3 = vf_2
    a0_3 = af2
    vf_3 = 0.0
    af3 = 0.0

    pitch_3, _, _, _ = get_xz_w_initial_conds(p2, p3, p0, p0, t0, tf3, v0_3, vf_3, a0_3, af3, dt)

    return pitch_1, pitch_2, pitch_3


def compute_quaternion_traj(roll_traj, euler_angle):
    # input euler angle is in X, Y, Z
    quaternion = np.zeros((roll_traj.shape[0], 4))  # Should be X, Y , Z, W

    # First fix the roll trajectory. In the euler angles there should be a sign change between
    # -180 -> -220, making it
    fixed_roll_traj = roll_traj.copy()
    for i in range(fixed_roll_traj.shape[0]):
        if fixed_roll_traj[i] < -180:
            # print(f"Found trajectory less than -180 i={i}, {fixed_roll_traj[i]} -> {360+fixed_roll_traj[i]}")
            fixed_roll_traj[i] = 360+fixed_roll_traj[i]

    for i in range(roll_traj.shape[0]):
        temp_euler = R.from_euler('xyz', [fixed_roll_traj[i], euler_angle[1], euler_angle[2]], degrees=True)
        quat = temp_euler.as_quat()
        quaternion[i, :] = quat

    return quaternion


def compute_x_z_angle(pos_start, pos_back, pos_front, pos_center, pos_last, angle_steps, euler_angle, dt):
    tf1 = 1.1
    tf2 = 1.5
    tf3 = 1.3
    tf4 = 0.5

    xnew, xnew2, xnew3, xnew4 = compute_x(pos_start, pos_back, pos_front, pos_center, pos_last, tf1, tf2, tf3, tf4, dt)

    znew, znew3, znew4 = compute_z(pos_start, pos_back, pos_front, pos_center, pos_last, tf1, tf3, tf4, dt)
    znew2 = xnew2.copy()
    znew2[:] = pos_front[0][2]

    pitch_1, pitch_2, pitch_3 = compute_angle(angle_steps, tf1, tf2, tf3, dt)
    pitch_4 = znew4.copy()
    pitch_4[:] = angle_steps[-1]

    x = np.concatenate((xnew[:-1], xnew2[:-1], xnew3[:-1], xnew4))
    z = np.concatenate((znew[:-1], znew2[:-1], znew3[:-1], znew4))
    roll_traj = np.concatenate((pitch_1[:-1], pitch_2[:-1], pitch_3[:-1], pitch_4))

    quaternion_traj = compute_quaternion_traj(roll_traj, euler_angle)

    # Create the quintic trajectory with the orientation
    quintic_traj = np.zeros((xnew.shape[0] + xnew2.shape[0] + xnew3.shape[0] + xnew4.shape[0] - 3, 7))  # For

    quintic_traj[:, 0] = x
    quintic_traj[:, 1] = pos_front[0][1]
    quintic_traj[:, 2] = z
    quintic_traj[:, 3] = quaternion_traj[:, 0]
    quintic_traj[:, 4] = quaternion_traj[:, 1]
    quintic_traj[:, 5] = quaternion_traj[:, 2]
    quintic_traj[:, 6] = quaternion_traj[:, 3]

    tf = tf1 + tf2 + tf3 + tf4

    return quintic_traj, tf


def generate_fling_trajectory(pose_left, pose_right, euler_init, dt, fling_height, grasp_height,
                              y_fixed=True, forward=True):
    """
    :param pose_left: Initial position of the left picker. np.array, shape(3,)
    :param pose_right: Initial position of the right picker. np.array, shape(3,)
    :param dt: dt of the simulator to compute the length of the trajectory points
    :param modify_x: whether to apply the trajectory in the X axis fixing the values of Y to the initial position.
    :return:
    """

    if y_fixed:
        x_start = pose_left[0]
        y_start = pose_left[1]
    else:  # We need to use X for perform the computation, we will change the axis back after the computation
        x_start = pose_left[1]
        y_start = pose_left[0]

    if forward:  # The robot moves on the indicated axis forward, else it moves first backwards
        fwd_sign = 1
    else:
        fwd_sign = -1

    z_start = pose_left[2]
    z_middle = fling_height
    z_end = grasp_height


    pos_start = [[x_start, y_start, z_start]]
    # pos_back = [[x_start-fwd_sign*0.15, y_start, z_middle]]
    pos_back = [[x_start - fwd_sign * 0.25, y_start, z_middle]]
    pos_front = [[x_start+fwd_sign*0.3, y_start, z_middle]]  # It should go at least up to 0.4 be 0.64
    pos_center = [[x_start+fwd_sign*0.15, y_start, z_end]]
    pos_last = [[x_start-fwd_sign*0.1, y_start, z_end]]



    # From -170 -> -145 -> 121 --> -170
    angle_steps = [euler_init[0], euler_init[0]+59, euler_init[0]-59, euler_init[0], euler_init[0]]

    quintic_traj, tf = compute_x_z_angle(pos_start, pos_back, pos_front, pos_center, pos_last,
                                         angle_steps, euler_init, dt)
    quintic_left = quintic_traj.copy()
    quintic_right = quintic_traj.copy()

    if y_fixed:
        quintic_right[:, 1] = pose_right[1]
    else:
        quintic_left[:, 0] = pose_left[0]
        quintic_left[:, 1] = quintic_traj[:, 0]
        quintic_right[:, 1] = quintic_traj[:, 0]
        quintic_right[:, 0] = pose_right[0]

    return quintic_left, quintic_right, tf

