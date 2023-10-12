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
   
   joint_ori = traj[0]

   #first apply joint movement if it is far from desired location (large motion not doable in single linear motion)
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin


   action = input("Move relative (Y/N)?").capitalize()
   delta = 5 #how many cm movement in x direction
   #franka.move_relative(params=[0.03, 0.00, 0.00], traj_duration=3.0*1.6) #for joint movement to origin
   franka.move_relative(params=[0.10, 0.00, 0.00], traj_duration=10.0*0.1) #for joint movement to origin

   # while action == 'Y':
   #    delta = 5 #how many cm movement in x direction
   #    franka.move_relative(params=[delta*0.01, 0.00, 0.00], traj_duration=delta) #for joint movement to origin
   #    time.sleep(delta)
   #    action = input("Move relative (Y/N)?").capitalize()