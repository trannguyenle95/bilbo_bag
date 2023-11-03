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
      width = 0.44
      V_max = 17.7
      A_max = 462
   elif args.Bag == "B":
      width = 0.37
      V_max = 17.7 #TODO: change
      A_max = 462 #TODO: change
   #TODO: add other bags (maybe used pandas df?)
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


   A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=True)
   print("Area (cm2): ", A_poly_rim, "Vol (l): ", Vol, " E_rim :", E_rim)

   data = {
   "A_CH_rim": [A_CH_rim],
   "A_poly_rim": [A_poly_rim],
   "Vol": [Vol],
   "E_rim": [E_rim],
   "Action": ["initial state"]
   }
   df = pd.DataFrame(data)


   input("Perform dynamic primitive")
   time.sleep(5.0) #add 5s sleep to manually move bag
   actions += 1
   #time.sleep(10) #sleep 10s when operating robots alone
   franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
   time.sleep(tf) #let motion finish before plotting and closing grippers
   real_pose = np.genfromtxt(pose_file, delimiter=',') #used to know the xdist at the end of a flip
   A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=True)
   print("Area (cm2): ", A_poly_rim, "Vol (l): ", Vol, " E_rim :", E_rim)
   new_pose = real_pose[-1]
   x_min = real_pose[-1][0]
   x_max = 0.65

   print("actions: ", actions)

   df.loc[len(df.index)] = [A_CH_rim, A_poly_rim, Vol, E_rim, "F"] 

   while actions <= max_actions:
      actions += 1
      if (Vol < 0.5*V_max) or (A_poly_rim < 0.5*A_max):
         action = "F"
         print("ACTION: ", action)
         franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
         time.sleep(3.0)
         franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
         time.sleep(tf)
         real_pose = np.genfromtxt(pose_file, delimiter=',')
         new_pose = real_pose[-1]

      elif E_rim < 0.8:
         delta = 0.01 #how many cm movement in x direction
         if new_pose[0] + delta < x_max:
            action = "DI"
            print("ACTION: ", action)
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            time.sleep(0.5) #NOTE: need higher sleep time if I want to test with remote bag!
            new_pose[0] = new_pose[0] + delta
         else:
            #TODO: <<< Fix this so robot doesnt go through actions quickly here!!
            print("max xdist reached")
            print("pose: ", new_pose)
   
      elif E_rim > 1.2:
         delta = -0.01 #how many cm movement in x direction
         if new_pose[0] - delta > x_min:
            action = "DD"
            print("ACTION: ", action)
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            time.sleep(0.5) #NOTE: need higher sleep time if I want to test with remote bag!
            new_pose[0] = new_pose[0] + delta
         else:
            #TODO: <<< Fix this so robot doesnt go through actions quickly here!!
            print("max xdist reached")
            print("pose: ", new_pose)

      else:
         print("sufficient bag state reached in ", actions, "actions")
         break
      
      A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=True)
      print("Area (cm2): ", A_poly_rim, "Vol (l): ", Vol, " E_rim :", E_rim)
      print("actions: ", actions)

      df.loc[len(df.index)] = [A_CH_rim, A_poly_rim,  Vol, E_rim, action] 
      time.sleep(5.0) #add 5s sleep to manually move bag
      
   print("final state:")
   A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=True)
   print("Vol (l): ", Vol, " E_rim :", E_rim)
   #NOTE: previous loop breaks when sufficient state is reached, so get final state here + SHOW PLOT

   path = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'runs', args.Bag+'_'+args.DMP+"_"+args.Demo[:-len('.csv')]+'_'+str(args.Run)+'.csv')
   df.to_csv(path)


   franka.open_grippers_middle()