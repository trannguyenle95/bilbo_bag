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
   franka = Franka(init_node=True)
   franka.rate.sleep()

   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')

   #Import traj and duration from CSV
   traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_2h_flip2.csv", delimiter=',') #NOTE: set name here!
   vel_traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_vel_2h_flip2.csv", delimiter=',') #NOTE: set name here!

   dt = 0.001 #1/30 #/ 120 #NOTE: modify this to match generated CSV or set this FPS when generating csv!
   #NOTE: seems like I have to use higher dt than 1/FPS, even when there is a sleep in franka.py

   input("Move robots to origin")
   
   joint_ori = traj[0]

   #first apply joint movement if it is far from desired location (large motion not doable in single linear motion)
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin


   franka.close_grippers_middle()
   input("Close grippers")
   franka.close_grippers() #NEW
   

   tf = traj.shape[0] * dt

   print('tf:', tf)

   filepath = os.path.join(datafolder+"/"+"executed_joint_trajectory.csv") #NOTE: new file for joint traj!

   if os.path.exists(filepath):
      os.remove(filepath) 


   pose_file = os.path.join(datafolder+"/"+"joint_control_pose.csv")
   if os.path.exists(pose_file):
      os.remove(pose_file) 


   joint_vel_file = os.path.join(datafolder+"/"+"executed_joint_velocities.csv")
   if os.path.exists(joint_vel_file):
      os.remove(joint_vel_file) 


   input("Perform dynamic primitive")


   franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)

   time.sleep(tf) #let motion finish before plotting and closing grippers

   real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!


   #Plot pose from executed actual motion
   ref_pose = np.genfromtxt(datafolder+"/trajectories/"+"pose_2h_flip2.csv", delimiter=',') #NOTE: set name here!
   real_pose = np.genfromtxt(pose_file, delimiter=',') #NOTE: set name here!
   real_vel = np.genfromtxt(joint_vel_file, delimiter=',') #NOTE: set name here!
   
   real_acc = np.diff(real_vel) / 0.001; 
   print("real acc shape:", real_acc.shape)

   #TODO: add more plots
   plt.figure(1)
   plt.plot(traj[:, 0], '-', label="j1")
   plt.plot(traj[:, 1],'-', label="j2")
   plt.plot(traj[:, 2], '-', label="j3")
   plt.plot(traj[:, 3], '-', label="j4")
   plt.plot(traj[:, 4], '-', label="j5")
   plt.plot(traj[:, 5], '-', label="j6")
   plt.plot(traj[:, 6], '-', label="j7")
   plt.title("Position components")
   plt.plot(real_traj[:, 0], '--', label="j1")
   plt.plot(real_traj[:, 1], '--', label="j2")
   plt.plot(real_traj[:, 2], '--', label="j3")
   plt.plot(real_traj[:, 3], '--', label="j4")
   plt.plot(real_traj[:, 4], '--', label="j5")
   plt.plot(real_traj[:, 5], '--', label="j6")
   plt.plot(real_traj[:, 6], '--', label="j7")
   plt.legend()
   plt.show()

   #PLOT Cartesian Pose during joint control 
   plt.figure(2)
   plt.plot(ref_pose[:, 0], 'r-', label="ref pos_x")
   plt.plot(ref_pose[:, 1],'g-', label="ref pos_y")
   plt.plot(ref_pose[:, 2], 'b-', label="ref pos_z")
   plt.title("Position components")
   plt.plot(real_pose[:, 0], 'r--', label="actual pos_x")
   plt.plot(real_pose[:, 1], 'g--', label="actual pos_y")
   plt.plot(real_pose[:, 2], 'b--', label="actual pos_z")
   plt.legend()

   plt.figure(3)
   plt.plot(ref_pose[:, 3], 'r-', label="ref qx")
   plt.plot(ref_pose[:, 4], 'g-', label="ref qy")
   plt.plot(ref_pose[:, 5], 'b-', label="ref qz")
   plt.plot(ref_pose[:, 6], 'm-', label="ref qw")
   plt.title("Quaternion components")
   plt.plot(real_pose[:, 3], 'r--', label="actual qx")
   plt.plot(real_pose[:, 4], 'g--', label="actual qy")
   plt.plot(real_pose[:, 5], 'b--', label="actual qz")
   plt.plot(real_pose[:, 6], 'm--', label="actual qw")
   plt.legend()

   plt.figure(4)
   plt.title("Joint vels")
   plt.plot(real_vel[:, 0], '-', label="vel j1")
   plt.plot(real_vel[:, 1],'-', label="vel j2")
   plt.plot(real_vel[:, 2], '-', label="vel j3")
   plt.plot(real_vel[:, 3], '-', label="vel j4")
   plt.plot(real_vel[:, 4], '-', label="vel j5")
   plt.plot(real_vel[:, 5], '-', label="vel j6")
   plt.plot(real_vel[:, 6], '-', label="vel j7")
   plt.legend()
   plt.show()

   plt.figure(5)
   plt.title("Joint accs")
   plt.plot(real_acc[:, 0], '-', label="acc j1")
   plt.plot(real_acc[:, 1],'-', label="acc j2")
   plt.plot(real_acc[:, 2], '-', label="acc j3")
   plt.plot(real_acc[:, 3], '-', label="acc j4")
   plt.plot(real_acc[:, 4], '-', label="acc j5")
   plt.plot(real_acc[:, 5], '-', label="acc j6")
   plt.plot(real_acc[:, 6], '-', label="acc j7")
   plt.legend()
   plt.show()



   plt.show()

   franka.open_grippers()
   #franka.release_grippers() #NEW
   #franka.move(move_type='d',params=traj, traj_duration=4.0)
