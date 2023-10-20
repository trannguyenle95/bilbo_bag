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

sys.path.append(os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src'))

import SupportScripts.ROS_BagMetrics as BagMetrics

#matplotlib.use('TkAgg')

VELOCITY_TIME_MULTIPLIER = 1.0


if __name__ == '__main__':
   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')
   #create robot error flag files so that pipeline will not run if you forget to run robot stop script
   filepath = os.path.join(datafolder+"/"+"fr2_error.txt") #read errors from other robot
   f = open(filepath, 'x')
   f.close
   filepath = os.path.join(datafolder+"/"+"fr3_error.txt") #read errors from other robot
   f = open(filepath, 'x')
   f.close


   parser = argparse.ArgumentParser()
   parser.add_argument('Bag', type=str, help='Select bag: A-E')
   args = parser.parse_args()

   if args.Bag == "A":
      width = 0.44
      V_max = 17.7
      A_max = 462
   #TODO: add other bags (maybe used pandas df?)
   max_actions = 10
   actions = 0

   rospy.init_node('franka3_talker', anonymous=True)

   franka = Franka(topic='/franka/', node_name='franka2_3_talker')

   franka.rate.sleep()

   #Import traj and duration from CSV
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


   A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=False)
   print("A_CH_rim (cm2): ", A_CH_rim, "A_poly_rim (cm2): ", A_poly_rim, " Vol (l): ", Vol, " E_rim :", E_rim)

   input("Perform dynamic primitive")
   #time.sleep(10) #sleep 10s when operating robots alone

   franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)

   time.sleep(tf) #let motion finish before plotting and closing grippers

   real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!

   real_pose = np.genfromtxt(pose_file, delimiter=',') #NOTE: set name here!

   A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=False)
   print("A_CH_rim (cm2): ", A_CH_rim, "A_poly_rim (cm2): ", A_poly_rim, " Vol (l): ", Vol, " E_rim :", E_rim)

   new_pose = real_pose[-1]
   x_min = real_pose[-1][0]
   x_max = 0.65

   actions += 1
   print("actions: ", actions)

   while actions <= max_actions:
      if Vol < 0.5*V_max:
         print("redoing flip")
         joint_ori = traj[0]
         franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
         time.sleep(3.0)
         franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
         time.sleep(tf)
         real_pose = np.genfromtxt(pose_file, delimiter=',')
         new_pose = real_pose[-1]

      elif E_rim < 0.8:
         delta = 0.02 #how many cm movement in x direction
         if new_pose[0] + delta < x_max:
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            new_pose[0] = new_pose[0] + delta
         else:
            print("max xdist reached")
            print("pose: ", new_pose)
   
      elif E_rim > 1.2:
         delta = -0.02 #how many cm movement in x direction
         if new_pose[0] - delta > x_min:
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            new_pose[0] = new_pose[0] + delta
         else:
            print("max xdist reached")
            print("pose: ", new_pose)

      else:
         print("sufficient bag state reached in ", actions, "actions")
         break
      
      A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=False)
      print("A_CH_rim (cm2): ", A_CH_rim, "A_poly_rim (cm2): ", A_poly_rim, " Vol (l): ", Vol, " E_rim :", E_rim)

      actions += 1
      print("actions: ", actions)

   print("final state:")
   A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=True)
   print("A_CH_rim (cm2): ", A_CH_rim, "A_poly_rim (cm2): ", A_poly_rim, " Vol (l): ", Vol, " E_rim :", E_rim)

   franka.open_grippers_middle()