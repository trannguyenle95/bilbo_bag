import sys
sys.path.append("/home/david/catkin_ws/src/franka_david/scripts/qdc-manip")

import rospy
import time
import numpy as np
import cv2
from scipy import optimize
from copy import deepcopy
from scipy.spatial.transform import Rotation as R

from franka import (Franka, DEFAULT_ORN)
from envs.primitives import generate_fling_trajectory

import matplotlib
from matplotlib import pyplot as plt

matplotlib.use('TkAgg')

VELOCITY_TIME_MULTIPLIER = 1.0


if __name__ == '__main__':
    # recorder = Recorder()
    # Current set-up moves it in Y
    # New minimum Z = -0.184, Start Z=0.51419
    # Start pose = 0.35415; -0.0012209; 0.51051
    # End pose = 0.35415; -0.0042028; -0.20033
    rospy.init_node('franka3_talker', anonymous=True)

    joint_ori = [0.10658619127984564, -0.1326096948608778, -0.13426805117803325, -1.7488134096131844,
                  0.01546355958738321, 1.6926563864284092, 0.842766618054774]

    franka = Franka(topic='/franka/', node_name='franka2_3_talker')
    franka.rate.sleep()



    # Orientation for place on table = 0.99981; -0.01426; 0.010089; 0.0089319, might need to change Z
    # Initial orientation = 0.99834; -0.037063; 0.038108; -0.021968
    # Robot left forward is = 0.39081; -0.25195; 0.30031
    # Robot backward is 0.50887; 0.29526; 0.29362
    input("Move robots to origin")
    franka.move(move_type='j', params=joint_ori, traj_duration=3.0)

    pose_left = [0.53, -0.0012209, 0.51051]
    pose_right = [0.53, -0.0012209, 0.51051]

    # Defined as X, Y, Z
    euler_init = [-177.31653135,   -4.27029547,   -4.35228107]  # Before = [-179.24297822, -2.93434324, 42.65659799]
    # euler_init = [-177.31653135, -177.31653135, -177.31653135]

    fling_height = pose_left[2] + 0.06  # Chequered rag
    # fling_height = pose_left[2] + 0.16  # Maximum?
    # fling_height = pose_left[2] - 0.1  # Does not unfold
    grasp_height = 0.0  # Should be -0.20
    # grasp_height = 0.0  # Should be -0.20
    dt = 1 / 1000

    quintic_left, quintic_right, tf = generate_fling_trajectory(pose_left, pose_right, euler_init, dt,
                                                                fling_height, grasp_height,
                                                                y_fixed=False, forward=False)

    quintic_traj = quintic_left.copy()
    # Let's test first the trajectory with fixed orientation
    # quintic_traj[:, 0] = quintic_traj[0, 0]
    # quintic_traj[:, 1] = quintic_traj[0, 1]
    # quintic_traj[:, 3] = quintic_traj[0, 3]
    # quintic_traj[:, 4] = quintic_traj[0, 4]
    # quintic_traj[:, 5] = quintic_traj[0, 5]
    # quintic_traj[:, 6] = quintic_traj[0, 6]
    print(f"Quaternion x={quintic_traj[0, 3]}, y={quintic_traj[0, 4]}, z={quintic_traj[0, 5]}, w={quintic_traj[0, 6]}")

    # Back orientation is = 0.95355; -0.037868; 0.049529; -0.2947
    # Forward orientation is = 0.87172; -0.052444; 0.0078392; 0.48713

    input("Perform dynamic primitive")
    franka.movedynamic_ori(quintic_traj=quintic_traj, tf=tf)


