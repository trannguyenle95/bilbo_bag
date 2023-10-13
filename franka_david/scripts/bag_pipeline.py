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

import os

matplotlib.use('TkAgg')

VELOCITY_TIME_MULTIPLIER = 1.0


if __name__ == '__main__':
   rospy.init_node('franka3_talker', anonymous=True)

   franka = Franka(topic='/franka/', node_name='franka2_3_talker')

   franka.rate.sleep()

   #Import traj and duration from CSV
   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')
   traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_demoDMP.csv", delimiter=',') #NOTE: set name here!
   vel_traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_vel_demoDMP.csv", delimiter=',') #NOTE: set name here!

   dt = 0.001

   action = input("Move relative (Y/N)?").capitalize()
   delta = 5 #how many cm movement in x direction
   #franka.move_relative(params=[0.03, 0.00, 0.00], traj_duration=3.0*1.6) #for joint movement to origin
   franka.move_relative(params=[0.10, 0.00, 0.00], traj_duration=10.0*0.1) #for joint movement to origin

   # input("Move robots to origin")

   # joint_ori = traj[0]
   # franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin


   # franka.close_grippers_middle()
   # input("Close grippers")
   # franka.close_grippers()

   # tf = traj.shape[0] * dt
   # print('tf:', tf)


   # filepath = os.path.join(datafolder+"/"+'executed_joint_trajectory.csv')
   # if os.path.exists(filepath):
   #    os.remove(filepath)

   # pose_file = os.path.join(datafolder+"/"+"joint_control_pose.csv")
   # if os.path.exists(pose_file):
   #    os.remove(pose_file) 

   # input("Perform dynamic primitive")
   # time.sleep(10) #sleep 10s when operating robots alone

   # franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)

   # time.sleep(tf) #let motion finish before plotting and closing grippers

   # real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!

   # real_pose = np.genfromtxt(pose_file, delimiter=',') #NOTE: set name here!

   # init_pose = real_pose[-1]

   # print("init pose: ", init_pose)

   # x_dist = init_pose[0]

   # action = input("Repeat (R), adjust distance (D), or stop (any other key)?").capitalize()

   # while action == 'R' or action == 'D':
   #    if action == 'R':
   #       joint_ori = traj[0]
   #       franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
   #       time.sleep(3.0)
   #       franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
   #       time.sleep(tf)
   #       real_pose = np.genfromtxt(pose_file, delimiter=',')
   #       init_pose = real_pose[-1]
   #       x_dist = init_pose[0]
   #    elif action == 'D':
   #       delta = 5 #how many cm movement in x direction
   #       if x_dist < 0.6:
   #          #franka.move_relative(params=[delta*0.01, 0.00, 0.00], traj_duration=5) #for joint movement to origin
   #          #time.sleep(1)
   #          franka.move_relative(params=[0.10, 0.00, 0.00], traj_duration=10.0) #for joint movement to origin
   #       else:
   #          print("max xdist reached")
   #          print("pose: ", pose)

   #    action = input("Repeat (R), adjust distance (D), or stop (any other key)?").capitalize()



   # franka.open_grippers_middle()