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


def distance_increase(refinement_actions, franka, x_curr, delta):
   refinement_actions += 1
   print("action: DI")
   action = "DI"
   franka.move_relative(params=[delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
   time.sleep(0.5) #NOTE: need higher sleep time if I want to test with remote bag!
   x_curr = x_curr + delta
   return action, refinement_actions, x_curr

def distance_decrease(refinement_actions, franka, x_curr, delta):
   refinement_actions += 1
   print("action: DD")
   action = "DD"
   franka.move_relative(params=[-delta, 0.00, 0.00], traj_duration=0.5) #for joint movement to origin
   time.sleep(0.5) #NOTE: need higher sleep time if I want to test with remote bag!
   x_curr = x_curr - delta
   return action, refinement_actions, x_curr

if __name__ == '__main__':
   print(">>>>> REMEMBER TO RUN ROBOT STOPPER SCRIPTS <<<<<")

   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')

   parser = argparse.ArgumentParser()
   parser.add_argument('Bag', type=str, help='Select bag: A-E')
   parser.add_argument('InitialState', type=str, help='Select whether initial state is "Easy" (not crumpled) or "Hard" (crumpled)?')
   parser.add_argument('Run', type=int, help='Write index of run with this bag/dmp/demo combination')
   args = parser.parse_args()

   if args.Bag == "A":
      V_max = 6.0
      A_max = 220
      init_dist = 0.573 #0.573 measured with franka_vizualization
   elif args.Bag == "B":
      V_max = 7.5
      A_max = 360
      init_dist = 0.558 # 0.558 measured with franka_vizualization
   elif args.Bag == "C":
      V_max = 12.5
      A_max = 370
      init_dist = 0.553 #0.553 measured with franka_vizualization
   elif args.Bag == "D":
      V_max = 14.0
      A_max = 530
      init_dist = 0.5189 #0.5189 measured with franka_vizualization
   elif args.Bag == "E":
      V_max = 22.0
      A_max = 550
      init_dist  = 0.499 #0.499 measured with franka_vizualization
   max_actions = 20
   refinement_actions = 0

   #Import traj and duration from CSV
   traj = np.genfromtxt(datafolder+"/trajectories/"+args.Bag+"_"+"Opt_DMP"+"_"+"joint_"+"10l_bag_flip.csv", delimiter=',')  #use this just for moving to origin!
  
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


   input("Start adjustment")
   delta = 0.01 #how many cm movement in x direction
   x_curr = init_dist
   x_min = init_dist
   x_max = 0.68


   while refinement_actions < max_actions:
      print("current xdist: ", x_curr)
      if ((A_alpha_rim < 0.6*A_max) or (Vol < 0.7*V_max)):

         if (E_rim) >= 0 and (E_rim < 1.0):
            if x_curr + delta < x_max:
               action, refinement_actions, x_curr = distance_increase(refinement_actions, franka, x_curr, delta)
            else:
               print("MAX xdist reached - decreasing dist instead!")
               action, refinement_actions, x_curr = distance_decrease(refinement_actions, franka, x_curr, delta)
      
         else:
            if x_curr - delta > x_min:
               print("action: DD")
               action, refinement_actions, x_curr = distance_decrease(refinement_actions, franka, x_curr, delta)
            else:
               print("MIN xdist reached - increasing dist instead!")
               action, refinement_actions, x_curr = distance_increase(refinement_actions, franka, x_curr, delta)

      else:
         print("Sufficient bag state reached in ", refinement_actions, "actions. Final state: ")
         break
      
      A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=False)
      print("Area %: ", A_alpha_rim/A_max, "Vol %: ", Vol/V_max, " E_rim:", E_rim)
      print("actions: ", refinement_actions)

      df.loc[len(df.index)] = [A_CH_rim, A_alpha_rim,  Vol, E_rim, action] 

   A_CH_rim, A_alpha_rim, Vol, E_rim = BagMetrics.calculate_metrics(args.Bag, displayPlot=True)
   #NOTE: previous loop breaks when sufficient state is reached, so get final state here + SHOW PLOT

   path = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'runs', 'refine', args.Bag+'_'+args.InitialState+str(args.Run)+'.csv')
   df.to_csv(path)


   open_grippers_msg = input("Open grippers (Y/N)?").upper()
   if(open_grippers_msg == "Y"):
      franka.release_grippers()