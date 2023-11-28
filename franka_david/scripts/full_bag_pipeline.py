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


def distance_increase(refinement_actions, franka, new_pose, delta):
   refinement_actions += 1
   print("action: DI")
   action = "DI"
   franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
   time.sleep(0.5) #NOTE: need higher sleep time if I want to test with remote bag!
   new_pose[0] = new_pose[0] + delta
   return action, refinement_actions, new_pose


def distance_decrease(refinement_actions, franka, new_pose, delta):
   if new_pose[0] - delta > x_min:
      refinement_actions += 1
      print("action: DD")
      action = "DD"
      franka.move_relative(params=[-delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
      time.sleep(0.5) #NOTE: need higher sleep time if I want to test with remote bag!
      new_pose[0] = new_pose[0] - delta
      return action, refinement_actions, new_pose


if __name__ == '__main__':
   print(">>>>> REMEMBER TO RUN ROBOT STOPPER SCRIPTS <<<<<")

   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')

   parser = argparse.ArgumentParser()
   parser.add_argument('Bag', type=str, help='Select bag: A-E')
   parser.add_argument('DMP', type=str, help='Select DMP version: tau_DMP / TC_DMP / Opt_DMP')
   parser.add_argument('Demo', type=str, help='Write name of demo file')
   parser.add_argument('InitialState', type=str, help='Select whether initial state is "Easy" (not crumpled) or "Hard" (crumpled)?')
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
   max_flips = 10
   max_refinement_actions = 20
   flip_actions = 0
   refinement_actions = 0

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
   flip_actions += 1
   franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
   time.sleep(tf) #let motion finish before plotting and closing grippers
   action = "F"
   real_pose = np.genfromtxt(pose_file, delimiter=',') #used to know the xdist at the end of a flip
   A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=False)
   print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)
   new_pose = real_pose[-1]
   x_min = real_pose[-1][0]
   x_max = 0.68

   print("flip actions: ", flip_actions)
   action = "F"
   df.loc[len(df.index)] = [A_CH_rim, A_alpha_rim,  Vol, E_rim, action] 

   delta = 0.01 #how many cm movement in x direction

   while flip_actions < max_flips:
      terminate = input("Did bag flip over (Y/N)?").upper()
      if (terminate == "Y"):
         path = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'runs', 'full_pipeline',args.DMP, args.Bag, args.InitialState, args.Bag+'_'+args.DMP+"_"+args.Demo[:-len('.csv')]+'_'+args.InitialState+str(args.Run)+'.csv')
         df.to_csv(path)
         sys.exit(0)

      if ((A_alpha_rim < 0.6*A_max) or (Vol < 0.7*V_max)):
         flip_actions += 1
         print("action: F")
         action = "F"
         joint_ori = traj[0]
         franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
         time.sleep(3.0)
         franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)
         time.sleep(tf)
         real_pose = np.genfromtxt(pose_file, delimiter=',')
         new_pose = real_pose[-1]

      else:
         break
      
      A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=False)
      print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)
      print("flip actions: ", flip_actions)
      df.loc[len(df.index)] = [A_CH_rim, A_alpha_rim,  Vol, E_rim, action] 

   while refinement_actions < max_refinement_actions:
      if ((A_alpha_rim >= 0.6*A_max) and (Vol >= 0.7*V_max) and (abs(E_rim) < 0.8) and (E_rim >= 0)):
         if new_pose[0] + delta < x_max:
            action, refinement_actions, new_pose = distance_increase(refinement_actions, franka, new_pose, delta)
         else:
            print("MAX xdist reached")
            break

   
      elif ((A_alpha_rim >= 0.6*A_max) and (Vol >= 0.7*V_max) and (abs(E_rim) < 0.8) and (E_rim < 0)):
         if new_pose[0] - delta > x_min:
               action, refinement_actions, new_pose = distance_decrease(refinement_actions, franka, new_pose, delta)
         else:
            print("MIN xdist reached")
            break

      elif ((A_alpha_rim < 0.6*A_max) or (Vol < 0.7*V_max)):
         if action == "DI":
            if new_pose[0] - 2*delta > x_min: #try to step back by 2*delta for larger margin to possible bad state
               print("undoing previous refinement action with margin")
               action, refinement_actions, new_pose = distance_decrease(refinement_actions, franka, new_pose, 2*delta)
            else:
               print("undoing previous refinement action")
               action, refinement_actions, new_pose = distance_decrease(refinement_actions, franka, new_pose, delta)
         else:
            if new_pose[0] + 2*delta < x_max: #try to step back by 2*delta for larger margin to possible bad state
               print("undoing previous refinement action with margin")
               action, refinement_actions, new_pose = distance_increase(refinement_actions, franka, new_pose, 2*delta)
            else:
               print("undoing previous refinement action")
               action, refinement_actions, new_pose = distance_increase(refinement_actions, franka, new_pose, delta)
         if refinement_actions == (max_refinement_actions - 1):
            #terminate run here so it is not risked that the final action moves back into a bad state if this one fixes it
            print("terminating so that the last action does not risk taking the bag to a worse state")
            A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=False)
            print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)
            print("refinement actions: ", refinement_actions)
            df.loc[len(df.index)] = [A_CH_rim, A_alpha_rim,  Vol, E_rim, action] 
            break
      else:
         print("Sufficient bag state reached in ", refinement_actions, "refinement actions. Final state: ")
         break
      
      A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=False)
      print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)
      print("refinement actions: ", refinement_actions)
      df.loc[len(df.index)] = [A_CH_rim, A_alpha_rim,  Vol, E_rim, action] 

   A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=True)
   #NOTE: previous loop breaks when sufficient state is reached, so get final state here + SHOW PLOT

   path = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'runs', 'full_pipeline',args.DMP, args.Bag, args.InitialState, args.Bag+'_'+args.DMP+"_"+args.Demo[:-len('.csv')]+'_'+args.InitialState+str(args.Run)+'.csv')
   df.to_csv(path)





   open_grippers_msg = input("Open grippers (Y/N)?").upper()
   if(open_grippers_msg == "Y"):
      franka.release_grippers()