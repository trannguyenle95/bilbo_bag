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

import pandas as pd 

#matplotlib.use('TkAgg')

VELOCITY_TIME_MULTIPLIER = 1.0


if __name__ == '__main__':
   print(">>>>> REMEMBER TO RUN ROBOT STOPPER SCRIPTS <<<<<")

   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')

   parser = argparse.ArgumentParser()
   parser.add_argument('Bag', type=str, help='Select bag: A-E')
   parser.add_argument('DMP', type=str, help='Select DMP version: tau_DMP / TC_DMP / Opt_DMP')
   parser.add_argument('Demo', type=str, help='Write name of demo file')
   parser.add_argument('Run', type=int, help='Write index of run with this bag/dmp/demo combination')
   args = parser.parse_args()

   if args.Bag == "A":
      V_max = 6.2
      A_max = 220
   elif args.Bag == "B":
      V_max = 7.5
      A_max = 380
   elif args.Bag == "C":
      V_max = 16.0
      A_max = 420
   elif args.Bag == "D":
      V_max = 16.0
      A_max = 590
   elif args.Bag == "E":
      V_max = 23.5
      A_max = 420
   max_actions = 20
   actions = 0

   #Import traj and duration from CSV
   traj = np.genfromtxt(datafolder+"/trajectories/"+args.Bag+"_"+args.DMP+"_"+"joint_"+args.Demo, delimiter=',')
   vel_traj = np.genfromtxt(datafolder+"/trajectories/"+args.Bag+"_"+args.DMP+"_"+"joint_vel_"+args.Demo, delimiter=',')

   rospy.init_node('franka3_talker', anonymous=True)
   franka = Franka(topic='/franka/', node_name='franka2_3_talker')
   franka.rate.sleep()

   dt = 0.001

   input("Move robots to origin")

   joint_ori = traj[0]
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin


   franka.close_grippers_middle()
   input("Close grippers")
   franka.close_grippers()

   tf = traj.shape[0] * dt
   print('tf:', tf)


   #used to know the xdist at the end of a flip
   pose_file = os.path.join(datafolder+"/"+"joint_control_pose.csv")
   if os.path.exists(pose_file):
      os.remove(pose_file)


   A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=False)
   print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)

   data = {
   "A_CH_rim": [A_CH_rim],
   "A_alpha_rim": [A_alpha_rim],
   "Vol": [Vol],
   "E_rim": [E_rim],
   "Action": ["initial state"]
   }
   df = pd.DataFrame(data)


   input("Perform dynamic primitive")
   time.sleep(5.0) #add 5.0 seconds wait to move the bag
   actions += 1
   franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
   time.sleep(tf) #let motion finish before plotting and closing grippers
   action = "F"
   real_pose = np.genfromtxt(pose_file, delimiter=',') #used to know the xdist at the end of a flip
   A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=False)
   print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)
   new_pose = real_pose[-1]
   x_min = real_pose[-1][0]
   x_max = 0.65

   print("actions: ", actions)
   action = "F"
   df.loc[len(df.index)] = [A_CH_rim, A_alpha_rim,  Vol, E_rim, action] 

   while actions < max_actions:
      actions += 1
      if ((A_alpha_rim < 0.5*A_max) or (Vol < 0.5*V_max)):
         print("action: F")
         time.sleep(5.0) #add 5.0 seconds wait to move the bag
         action = "F"
         joint_ori = traj[0]
         franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
         time.sleep(3.0)
         franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
         time.sleep(tf)
         real_pose = np.genfromtxt(pose_file, delimiter=',')
         new_pose = real_pose[-1]

      elif E_rim < 0.8:
         delta = 0.01 #how many cm movement in x direction
         if new_pose[0] + delta < x_max:
            print("action: DI")
            action = "DI"
            time.sleep(5.0) #add 5.0 seconds wait to move the bag
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            time.sleep(0.5) #NOTE: need higher sleep time if I want to test with remote bag!
            new_pose[0] = new_pose[0] + delta
         else:
            print("MAX xdist reached")
            print("pose: ", new_pose)
            time.sleep(5.0) #add 5.0 seconds wait to move the bag
   
      elif E_rim > 1.2:
         delta = -0.01 #how many cm movement in x direction
         if new_pose[0] - delta > x_min:
            print("action: DD")
            action = "DD"
            time.sleep(5.0) #add 5.0 seconds wait to move the bag
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            time.sleep(0.5) #NOTE: need higher sleep time if I want to test with remote bag!
            new_pose[0] = new_pose[0] + delta
         else:
            print("MIN xdist reached")
            print("pose: ", new_pose)
            time.sleep(5.0) #add 5.0 seconds wait to move the bag

      else:
         print("sufficient bag state reached in ", actions, "actions")
         break
      
      A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=False)
      print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)
      print("actions: ", actions)

      df.loc[len(df.index)] = [A_CH_rim, A_alpha_rim,  Vol, E_rim, action] 

   print("final state:")
   print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)
   A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=True)
   #NOTE: previous loop breaks when sufficient state is reached, so get final state here + SHOW PLOT

   path = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'runs', args.Bag+'_'+args.DMP+"_"+args.Demo[:-len('.csv')]+'_'+str(args.Run)+'.csv')
   df.to_csv(path)


   franka.open_grippers_middle()