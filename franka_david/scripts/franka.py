import numpy as np
import rospy

from time import sleep, time
from copy import deepcopy

from franka_david.msg import MovetoPy
from franka_david.msg import MotionPy
from franka_david.msg import JointMotionPy
from franka_david.msg import JointTrajPy #ADDED 29.08
from franka_david.msg import GripperPy
from franka_david.msg import GripperGraspPy
from franka_david.msg import RemoteMotionPy
from geometry_msgs.msg import PoseStamped


from primitives import get_xyz_coordinates_given_mid_theta_ref_real

DEFAULT_ORN = [0.0, np.pi / 4, 0.0, 0.0]  # Quaternion
# OUT_OF_WAY = [0.08788468, 0.576599, 0.2247749]  # With modified X=Z
OUT_OF_WAY_MIDDLE = [0.50826, 0.059283, 0.45]
OUT_OF_WAY_END = [0.03718, 0.50722, 0.45]

# JOINTS_CENTER = [0.0, 0.0, 0.0, -np.pi/2,
#                  0.0, 1.5274167884720695, 0.8371567886289623]
JOINTS_CENTER = [0.0, 0.0, 0.0, -np.pi/2,
                 0.0, 1.6244422215865644, 0.8171567886289623]
# This seems to be more centered
# (-0.06806987871617443, -0.048338811235469674, -1.6072872920287282,
#  -1.6072872920287282, -0.12901891718086764, 1.6244422215865644, -0.6810651746177011]

JOINTS_OUT_OF_WAY = [-1.5402842149198652, 0.0, 0.0, -np.pi/2,
                     0.0, 1.6244422215865644, 0.8171567886289623]

# WS_PC = [10, 370, 416, 675] # Old one
# WS_PC = [10, 370, 492, 831]
# WS PC MASK FOR 1080/1920
# WS_PC_MASK = [136, 715, 708, 1290]  # When rotated
# WS_PC = [201, 671, 858, 1161]  # first are limits in image Y

# WS MASKS FOR 720/1280
WS_PC_MASK = [91, 477, 472, 860]  # When rotated
WS_PC = [134, 447, 572, 774]  # first are limits in image Y

# Manipulation limits
# WS_PC = [65, 535, 150, 453]  # first are limits in image Y
FRANKA_TIME = 1

WORKSPACE_SURFACE_MIN = 0.105  # Maximum 0.764 with 140 // maximum 0.778 with 85 - 67 // 63 with cardboard
WORKSPACE_SURFACE_MAX = 0.5

GRIPPER_OPEN = 0.07
GRIPPER_CLOSED = 0.0006
GRIPPER_MIDDLE_OPEN = 0.05
GRIPPER_MIDDLE_CLOSED = 0.002 #originally 0.025
GRIPPER_SPEED = 0.05
GRIPPER_FORCE = 40  # Newtons

bottom_right = [0.72368, -0.35026, 0.12776]  # 78265126
top_left = [0.30053, 0.35026, 0.12776]
bottom_left = [0.30053, -0.35026, 0.12776]
top_right = [0.72368, 0.35026, 0.12776]  # 7761067


def bound_grasp_pos(pos, z_offset=0.08):
    # TODO Modify so that this is within the correct workspace limits
    pos = deepcopy(pos)
    # bound the X and Y
    pos[0] = max(top_left[0], pos[0])  # max 0.11 vs e.g. 0.08
    pos[0] = min(top_right[0], pos[0])  # min 0.87 vs e.g. 0.9
    pos[1] = max(bottom_left[1], pos[1])
    pos[1] = min(top_left[1], pos[1])  # min 0.45 vs e.g. 0.6
    # grasp slightly lower than detected depth
    pos[2] -= z_offset
    pos[2] = max(WORKSPACE_SURFACE_MIN, pos[2])
    pos[2] = min(WORKSPACE_SURFACE_MAX, pos[2])
    return pos


class Franka:

    def __init__(self, topic='/frankapy/', init_node=False, node_name='franka3'):
        self.ros_moveto_pub = rospy.Publisher(topic+'moveto', MovetoPy,
                                              queue_size=10)
        self.ros_motion_pub = rospy.Publisher(topic+'motion', MotionPy,
                                              queue_size=10)
        self.ros_remotemotion_pub = rospy.Publisher(topic+'remote_motion', RemoteMotionPy,
                                                    queue_size=10)
        self.ros_motion_ori_pub = rospy.Publisher(topic + 'motion_ori', MotionPy,
                                                    queue_size=10)
        self.ros_joint_pub = rospy.Publisher(topic+'jointmotion', JointMotionPy,
                                             queue_size=10)
        self.ros_jointTraj_pub = rospy.Publisher(topic+'joint_trajectory', JointTrajPy,
                                             queue_size=10) #ADDED 29.08
        self.ros_gripper_move_pub = rospy.Publisher(topic+'gripper_move', GripperPy,
                                                    queue_size=10)
        self.ros_gripper_grasp_pub = rospy.Publisher(topic+'gripper_grasp', GripperGraspPy,
                                                     queue_size=10)

        if init_node:
            rospy.init_node(node_name, anonymous=True)

        self.rate = rospy.Rate(10)  # 10hz

    def send_home(self):
        self.move(move_type='j', params=JOINTS_CENTER, traj_duration=0.0)
        sleep(5)

    def out_of_the_way(self):
        self.move(move_type='j', params=JOINTS_OUT_OF_WAY, traj_duration=0.0)
        sleep(5)

    def open_grippers(self):
        self.release_grippers()
        print("Opening grippers")
        gripper_msg = GripperPy()
        # Width in meters
        gripper_msg.width = GRIPPER_OPEN
        # Speed in m/s
        gripper_msg.speed = GRIPPER_SPEED
        gripper_msg.enable = True
        self.ros_gripper_move_pub.publish(gripper_msg)
        self.rate.sleep()
        sleep(0.5)

    def close_grippers(self):
        gripper_msg = GripperGraspPy()
        # Width in meters
        gripper_msg.distance = GRIPPER_CLOSED
        # Speed in m/s
        gripper_msg.speed = GRIPPER_SPEED
        gripper_msg.force = GRIPPER_FORCE
        gripper_msg.release = False
        gripper_msg.enable = True
        self.rate.sleep()
        self.ros_gripper_grasp_pub.publish(gripper_msg)
        self.rate.sleep()
        sleep(1)  # Let's wait one second for the gripper to close

    def release_grippers(self):
        print("Releasing grippers")
        gripper_msg = GripperGraspPy()
        # Width in meters
        gripper_msg.distance = 0.0
        # Speed in m/s
        gripper_msg.speed = 0.0
        gripper_msg.force = 0.0
        gripper_msg.release = True
        gripper_msg.enable = True
        self.rate.sleep()
        self.ros_gripper_grasp_pub.publish(gripper_msg)
        self.rate.sleep()
        # sleep(1)

    def close_grippers_middle(self):
        print("Opening grippers")
        gripper_msg = GripperPy()
        # Width in meters
        gripper_msg.width = GRIPPER_MIDDLE_CLOSED
        # Speed in m/s
        gripper_msg.speed = GRIPPER_SPEED
        gripper_msg.enable = True
        self.ros_gripper_move_pub.publish(gripper_msg)
        self.rate.sleep()
        sleep(0.5)
        # self.release_grippers()  # Let's make sure we release before anything

    def open_grippers_middle(self):
        self.release_grippers()  # Let's make sure we release before anything
        print("Opening grippers")
        gripper_msg = GripperPy()
        # Width in meters
        gripper_msg.width = GRIPPER_MIDDLE_OPEN
        # Speed in m/s
        gripper_msg.speed = GRIPPER_SPEED
        gripper_msg.enable = True
        self.ros_gripper_move_pub.publish(gripper_msg)
        self.rate.sleep()
        sleep(0.5)

    def compute_dynamic(self, init_pos, end_pos, mid_theta_ref):
        start_pos = np.array(init_pos)
        target_pos = np.array(end_pos)  # Get the target position
        # init_pos = self.action_tool._get_pos()[0]  # Get the initial position

        start_pos[0, 2] = init_pos[0][1]
        start_pos[0, 1] = init_pos[0][2]
        target_pos[0, 2] = end_pos[0][1]
        target_pos[0, 1] = end_pos[0][2]

        # Compute the cubic/quintic polynomial
        quintic_traj, theta, coefs, tf = get_xyz_coordinates_given_mid_theta_ref_real(start_pos, target_pos,
                                                                                      mid_theta_ref=mid_theta_ref)
        # Initial and ending Y position should be the grasp height
        # quintic_traj[:, 1] += end_pos[0][1]

        # Need to swap Y with Z
        swapped_quintic = quintic_traj[:, [0, 2, 1]]  # For
        # Initial and ending Y position should be the grasp height
        swapped_quintic[:, 2] = init_pos[0][2]

        return swapped_quintic, tf


    def movedynamic_ori(self, quintic_traj, tf):

        self.move(move_type='o', params=quintic_traj, traj_duration=tf)

    def movedynamic(self, quintic_traj, tf):

        import matplotlib
        matplotlib.use('TkAgg')
        # from matplotlib import pyplot as plt
        # figure, ax = plt.subplots(3, 1)
        # x = np.arange(swapped_quintic.shape[0])
        # ax[0].plot(x, swapped_quintic[:, 0], label="X")
        # ax[1].plot(x, swapped_quintic[:, 1], label="Y")
        # ax[2].plot(x, swapped_quintic[:, 2], label="Z")
        # plt.show()

        from matplotlib import pyplot as plt

        dt = 1/1000
        # t = np.linspace(start=0, step=dt, end=tf+dt)
        start = 0.0
        stop = tf#+dt
        step = dt
        t = np.linspace(start, stop, int((stop - start) / step + 2))
        # t = np.arange(quintic_traj.shape[0])
        x = quintic_traj[:, 0]
        print(f"Shapes x={x.shape[0]}, t={t.shape[0]}")
        x_vel = np.gradient(x, dt)
        x_acc = np.gradient(x_vel, dt)
        print(f"Minimum vel is ={x_vel.min()}, max vel = {x_vel.max()}),"
              f" min acc={x_acc.min()}, max acc = {x_acc.max()}")
        z = quintic_traj[:, 2]
        z_vel = np.gradient(z, dt)
        z_acc = np.gradient(z_vel, dt)
        # figure, ax = plt.subplots(2, 3)
        # ax[0, 0].plot(t, x)
        # ax[0, 0].set_ylabel("X (m)")
        # ax[0, 1].plot(t, x_vel)
        # ax[0, 1].set_ylabel("Velocity (m/s)")
        # ax[0, 2].plot(t, x_acc)
        # ax[0, 2].set_ylabel("Acceleration (m/s2)")
        # ax[1, 0].plot(t, z)
        # ax[1, 0].set_ylabel("Z (m)")
        # ax[1, 1].plot(t, z_vel)
        # ax[1, 1].set_ylabel("Velocity (m/s)")
        # ax[1, 2].plot(t, z_acc)
        # ax[1, 2].set_ylabel("Acceleration (m/s2)")
        # plt.suptitle(f'Time={tf}, min vel is ={x_vel.min()}, max vel = {x_vel.max()}')
        # plt.show()

        # TODO Double check that XYZ is correct
        self.move(move_type='d', params=quintic_traj, traj_duration=tf)


    # Example code for new xz-w initial conds
#pos_start = [[0.3, 0.30, 0.50]]
# pos_back = [[-0.1, 0.30, 0.50]]
# pos_front = [[0.64, 0.30, 0.50]]

# x0 = pos_start[0][0]
# xf = pos_back[0][0]
# z0 = 0.5
# zf = 0.5

# t0=0
# v0=0
# a0=0
# vf=-0.5
# af=0.0
# tf=3.0

# xnew, znew, _, _ = get_xz_w_initial_conds(x0, xf, z0, zf, t0, tf, v0, vf, a0, af)
#x0 = pos_back[0][0]
# xf = pos_front[0][0]
# z0 = 0.5
# zf = 0.5

# t0=0.0
# v0=-0.5
# a0=0
# vf=0.0
# af=0.0
# tf2=3.0

# xnew2, znew2, _, _ = get_xz_w_initial_conds(x0, xf, z0, zf, t0, tf2, v0, vf, a0, af)

    def movel(self, params, traj_duration: float = 10.0, **kwargs):
        # sleep the time of execution

        # TODO Check that params is a single list and not multiple lists
        self.move(move_type='l', params=params, traj_duration=traj_duration)
        # sleep for the amount of time the trajectory is executed to let the robot finish + 0.5 s
        # so that the trajectory finishes
        # sleep(traj_duration+0.5)

    def move(self, move_type, params, traj_duration: float, **kwargs):

        if move_type == 'l':  # Linear move
            params = params[0]
            assert traj_duration != 0.0
            pose_msg = MovetoPy()
            pose_msg.duration = traj_duration
            pose_msg.pose.pose.position.x = params[0]
            pose_msg.pose.pose.position.y = params[1]
            pose_msg.pose.pose.position.z = params[2]
            pose_msg.pose.pose.orientation.x = params[3]
            pose_msg.pose.pose.orientation.y = params[4]
            pose_msg.pose.pose.orientation.z = params[5]
            pose_msg.pose.pose.orientation.w = params[6]
            # pose_msg.pose.pose.orientation.x = DEFAULT_ORN[0]
            # pose_msg.pose.pose.orientation.y = DEFAULT_ORN[1]
            # pose_msg.pose.pose.orientation.z = DEFAULT_ORN[2]
            # pose_msg.pose.pose.orientation.w = DEFAULT_ORN[3]
            pose_msg.enable = True
            pose_msg.pose.header.seq = 1
            pose_msg.pose.header.stamp = rospy.Time.now()
            pose_msg.pose.header.frame_id = "move"
            # print("Sending new Cartesian position")
            # rospy.sleep(0.01) #originally 0.5 - keep it low here so actual dt is close to what is set in script!
            self.ros_moveto_pub.publish(pose_msg)
            sleep(traj_duration)

        elif move_type == 'd':  # Dynamic move
            assert traj_duration != 0.0
            pose_msg = MotionPy()
            pose_msg.duration = traj_duration
            trajectory = []
            for i in range(params.shape[0]):
                traj_i_msg = PoseStamped()

                traj_i_msg.pose.position.x = params[i, 0]
                traj_i_msg.pose.position.y = params[i, 1]
                traj_i_msg.pose.position.z = params[i, 2]
                traj_i_msg.pose.orientation.x = DEFAULT_ORN[0]
                traj_i_msg.pose.orientation.y = DEFAULT_ORN[1]
                traj_i_msg.pose.orientation.z = DEFAULT_ORN[2]
                traj_i_msg.pose.orientation.w = DEFAULT_ORN[3]

                traj_i_msg.header.seq = 1
                traj_i_msg.header.stamp = rospy.Time.now()
                traj_i_msg.header.frame_id = "move"
                trajectory.append(traj_i_msg)
            pose_msg.trajectory = trajectory
            pose_msg.enable = True
            print("Sending dynamic trajectory")
            rospy.sleep(0.1)
            self.ros_motion_pub.publish(pose_msg)
            sleep(traj_duration)

        elif move_type == 'o':  # Move trajectory with orientation
            assert traj_duration != 0.0
            pose_msg = MotionPy()
            pose_msg.duration = traj_duration
            pose_msg.enable = True
            trajectory = []
            for i in range(params.shape[0]):
                traj_i_msg = PoseStamped()

                traj_i_msg.pose.position.x = params[i, 0]
                traj_i_msg.pose.position.y = params[i, 1]
                traj_i_msg.pose.position.z = params[i, 2]
                traj_i_msg.pose.orientation.x = params[i, 3]
                traj_i_msg.pose.orientation.y = params[i, 4]
                traj_i_msg.pose.orientation.z = params[i, 5]
                traj_i_msg.pose.orientation.w = params[i, 6]

                traj_i_msg.header.seq = 1
                traj_i_msg.header.stamp = rospy.Time.now()
                traj_i_msg.header.frame_id = "move"
                trajectory.append(traj_i_msg)
            pose_msg.trajectory = trajectory
            print("Sending dynamic trajectory")
            rospy.sleep(0.1)
            self.ros_motion_ori_pub.publish(pose_msg)
            sleep(traj_duration)

        elif move_type == 'jt': #ADDED 29.08
            assert traj_duration != 0.0
            pose_msg = JointTrajPy()
            pose_msg.duration = traj_duration
            pose_msg.enable = True
            for i in range(params.shape[0]):
                pose_msg.joint0.append(params[i,0])
                pose_msg.joint1.append(params[i,1])
                pose_msg.joint2.append(params[i,2])
                pose_msg.joint3.append(params[i,3])
                pose_msg.joint4.append(params[i,4])
                pose_msg.joint5.append(params[i,5])
                pose_msg.joint6.append(params[i,6])

            print("Sending new joint trajectory")
            rospy.sleep(0.5)
            self.ros_jointTraj_pub.publish(pose_msg)
            sleep(2)

        elif move_type == 'j':
            pose_msg = JointMotionPy()
            pose_msg.joint0 = params[0]
            pose_msg.joint1 = params[1]
            pose_msg.joint2 = params[2]
            pose_msg.joint3 = params[3]
            pose_msg.joint4 = params[4]
            pose_msg.joint5 = params[5]
            pose_msg.joint6 = params[6]
            pose_msg.enable = True
            print("Sending new joint position")
            rospy.sleep(0.5)
            self.ros_joint_pub.publish(pose_msg)
            sleep(2)

        if move_type == 'r':  # Linear move
            params = params[0]
            requested_time = params[0]
            print(f"Requested time = {requested_time}")
            assert traj_duration != 0.0
            pose_msg = RemoteMotionPy()
            pose_msg.duration = traj_duration
            pose_msg.enable = True
            pose_msg.pose.header.seq = 1
            pose_msg.pose.header.stamp = requested_time
            pose_msg.pose.header.frame_id = "move"
            print("Sending new Cartesian position")
            rospy.sleep(0.5)
            self.ros_moveto_pub.publish(pose_msg)
            sleep(traj_duration)

        self.rate.sleep()
