import sys
from tkinter import E
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
import argparse

#sys.path.append(os.path.dirname(os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')))
#import ROS_BagMetrics

sys.path.append(os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src'))

import SupportScripts.ROS_BagMetrics as BagMetrics

#matplotlib.use('TkAgg')

VELOCITY_TIME_MULTIPLIER = 1.0


if __name__ == '__main__':
 
   parser = argparse.ArgumentParser()
   parser.add_argument('max_area', type=float, help='Maximum area of bag rim when opening is fully open and circular.')
   parser.add_argument('width', type=float, help='Maximum width of the bag.')
   args = parser.parse_args()

   rospy.init_node('franka3_talker', anonymous=True)

   franka = Franka(topic='/franka/', node_name='franka2_3_talker')

   franka.rate.sleep()

   #Import traj and duration from CSV
   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')
   traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_demoDMP.csv", delimiter=',') #NOTE: set name here!
   vel_traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_vel_demoDMP.csv", delimiter=',') #NOTE: set name here!

   dt = 0.001

   input("Move robots to origin")

   joint_ori = traj[0]
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin


   franka.close_grippers_middle()
   input("Close grippers")
   franka.close_grippers()

   tf = traj.shape[0] * dt
   print('tf:', tf)


   filepath = os.path.join(datafolder+"/"+'executed_joint_trajectory.csv')
   if os.path.exists(filepath):
      os.remove(filepath)

   pose_file = os.path.join(datafolder+"/"+"joint_control_pose.csv")
   if os.path.exists(pose_file):
      os.remove(pose_file) 

   input("Perform dynamic primitive")
   time.sleep(10) #sleep 10s when operating robots alone

   franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)

   time.sleep(tf) #let motion finish before plotting and closing grippers

   real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!

   real_pose = np.genfromtxt(pose_file, delimiter=',') #NOTE: set name here!

   # A_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.max_area, args.width, displayPlot=False)
   # print("A_rim (cm2): ", A_rim, " Vol (l): ", Vol, " E_rim :", E_rim)

   new_pose = real_pose[-1]
   x_min = real_pose[-1][0]
   x_max = 0.65

   action = input("Repeat (R), increase distance (DI), decrease distance (DD), or stop (any other key)?").upper()

   print("action: ", action)

   while action == 'R' or action == 'DI' or action == 'DD':
      if action == 'R':
         joint_ori = traj[0]
         franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
         time.sleep(3.0)
         franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
         time.sleep(tf)
         real_pose = np.genfromtxt(pose_file, delimiter=',')
         new_pose = real_pose[-1]
      elif action == 'DI':
         delta = 0.02 #how many cm movement in x direction
         if new_pose[0] + delta < x_max:
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            new_pose[0] = new_pose[0] + delta
         else:
            print("max xdist reached")
            print("pose: ", new_pose)
      elif action == 'DD':
         delta = -0.02 #how many cm movement in x direction
         if new_pose[0] - delta > x_min:
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            new_pose[0] = new_pose[0] + delta
         else:
            print("min xdist reached")
            print("pose: ", new_pose)

      # A_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.max_area, args.width, displayPlot=False)
      # print("A_rim (cm2): ", A_rim, " Vol (l): ", Vol, " E_rim :", E_rim)
      action = input("Repeat (R), increase distance (DI), decrease distance (DD), or stop (any other key)?").upper()



   franka.open_grippers_middle()