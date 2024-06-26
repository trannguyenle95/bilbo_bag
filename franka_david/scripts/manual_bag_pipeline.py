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
   parser.add_argument('OptiTrack', type=str, help='Gather OptiTrack metrics (Y/N)')
   parser.add_argument('Bag', type=str, help='Select bag: A-E')
   parser.add_argument('DMP', type=str, help='Select DMP version: tau_DMP / TC_DMP / Opt_DMP')
   parser.add_argument('Demo', type=str, help='Write name of demo file')
   parser.add_argument('Run', type=int, help='Write index of run with this bag/dmp/demo combination')
   args = parser.parse_args()

   if args.Bag == "A":
      V_max = 6.0
      A_max = 220
   elif args.Bag == "B":
      V_max = 7.5
      A_max = 360
   elif args.Bag == "C":
      V_max = 12.5
      A_max = 370
   elif args.Bag == "D":
      V_max = 14.0
      A_max = 530
   elif args.Bag == "E":
      V_max = 22.0
      A_max = 550
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


   open_grippers_msg = input("Open grippers (Y/N)?").upper()
   if(open_grippers_msg == "Y"):
      franka.release_grippers()
      #franka.close_grippers_middle()
      input("Close grippers")
      franka.close_grippers()
      print("CLOSING GRIPPERS IN (10s FR2 / 20s FR3)")
      time.sleep(20)

   tf = traj.shape[0] * dt
   print('tf:', tf)


   #used to know the xdist at the end of a flip
   pose_file = os.path.join(datafolder+"/"+"joint_control_pose.csv")
   if os.path.exists(pose_file):
      os.remove(pose_file)

   if args.OptiTrack == 'Y':
      A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=False)
      print("A_CH_rim (cm2): ", A_CH_rim, "A_poly_rim (cm2): ", A_poly_rim, " Vol (l): ", Vol, " E_rim :", E_rim)


      data = {
      "A_CH_rim": [A_CH_rim],
      "A_poly_rim": [A_poly_rim],
      "Vol": [Vol],
      "E_rim": [E_rim],
      "Action": ["initial state"]
      }
      df = pd.DataFrame(data)


   #for plotting
   # filepath = os.path.join(datafolder+"/"+"executed_joint_trajectory.csv") #NOTE: new file for joint traj!
   # if os.path.exists(filepath):
   #    os.remove(filepath) 

   input("Perform dynamic primitive")
   actions += 1
   #time.sleep(10) #sleep 10s when operating robots alone
   franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
   time.sleep(tf) #let motion finish before plotting and closing grippers
   real_pose = np.genfromtxt(pose_file, delimiter=',') #used to know the xdist at the end of a flip

   if args.OptiTrack == 'Y':
      A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=False)
      print("A_CH_rim (cm2): ", A_CH_rim, "A_poly_rim (cm2): ", A_poly_rim, " Vol (l): ", Vol, " E_rim :", E_rim)
   new_pose = real_pose[-1]
   x_min = real_pose[-1][0]
   x_max = 0.68

   print("actions: ", actions)


   if args.OptiTrack == 'Y':
      df.loc[len(df.index)] = [A_CH_rim, A_poly_rim,  Vol, E_rim, action] 
 

   #plotting
   # real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!

   # plt.figure(1)
   # plt.plot(traj[:, 0], '-', label="j1")
   # plt.plot(traj[:, 1],'-', label="j2")
   # plt.plot(traj[:, 2], '-', label="j3")
   # plt.plot(traj[:, 3], '-', label="j4")
   # plt.plot(traj[:, 4], '-', label="j5")
   # plt.plot(traj[:, 5], '-', label="j6")
   # plt.plot(traj[:, 6], '-', label="j7")
   # plt.title("Position components")
   # plt.plot(real_traj[:, 0], '--', label="j1")
   # plt.plot(real_traj[:, 1], '--', label="j2")
   # plt.plot(real_traj[:, 2], '--', label="j3")
   # plt.plot(real_traj[:, 3], '--', label="j4")
   # plt.plot(real_traj[:, 4], '--', label="j5")
   # plt.plot(real_traj[:, 5], '--', label="j6")
   # plt.plot(real_traj[:, 6], '--', label="j7")
   # plt.legend()
   # plt.show()

   # input("continue")


   #action = input("Repeat flip (F), increase distance (DI), decrease distance (DD), or stop (any other key)?").upper()
   action = 'F'


   print("action: ", action)

   delta = 0.01 #how many cm movement in x direction
   while action == 'F' or action == 'DI' or action == 'DD':
      if action == 'F':
         actions += 1
         joint_ori = traj[0]
         franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
         time.sleep(3.0)
         franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
         time.sleep(tf)
         real_pose = np.genfromtxt(pose_file, delimiter=',')
         new_pose = real_pose[-1]
      elif action == 'DI':
         actions += 1
         if new_pose[0] + delta < x_max:
            franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            new_pose[0] = new_pose[0] + delta
         else:
            print("max xdist reached")
            print("pose: ", new_pose)
      elif action == 'DD':
         actions += 1
         if new_pose[0] - delta > x_min:
            franka.move_relative(params=[-delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
            new_pose[0] = new_pose[0] - delta
         else:
            print("min xdist reached")
            print("pose: ", new_pose)

      if args.OptiTrack == 'Y':
         A_CH_rim, A_poly_rim, Vol, E_rim = BagMetrics.calculate_metrics(width, displayPlot=False)
         print("A_CH_rim (cm2): ", A_CH_rim, "A_poly_rim (cm2): ", A_poly_rim, " Vol (l): ", Vol, " E_rim :", E_rim)
      print("actions: ", actions)

      if args.OptiTrack == 'Y':
         df.loc[len(df.index)] = [A_CH_rim, A_poly_rim,  Vol, E_rim, action] 
      #action = input("Repeat flip (F), increase distance (DI), decrease distance (DD), or stop (any other key)?").upper()
      action = 'F'

   if args.OptiTrack == 'Y':
      #path = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'runs', args.Bag+'_'+args.DMP+"_"+args.Demo[:-len('.csv')]+'_'+str(args.Run)+'.csv')
      #df.to_csv(path)
      print("writing to file disabled")

   #franka.open_grippers_middle()

   
   open_grippers_msg = input("Open grippers (Y/N)?").upper()
   if(open_grippers_msg == "Y"):
      franka.release_grippers()