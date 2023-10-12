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
   # recorder = Recorder()
   # Current set-up moves it in Y
   # New minimum Z = -0.184, Start Z=0.51419
   # Start pose = 0.35415; -0.0012209; 0.51051
   # End pose = 0.35415; -0.0042028; -0.20033

   #NOTE: need to comment away something here?
   rospy.init_node('franka3_talker', anonymous=True)
   
   #For getting joint states:
   #roslaunch franka_visualization franka_visualization.launch robot_ip:=172.16.0.2 load_gripper:=True
   #Then run:
   #roslaunch franka_visualization franka_visualization.launch robot_ip:=172.16.0.2 load_gripper:=True
   
   #joint_ori = [-0.06153707569372121, 0.23268072162435294, -0.003733379824253959, -2.120620626949313, -0.07440938119840552, 2.374850448676014, 0.851590066155449] #35cm distance
   #joint_ori = [-0.04978088093284428, 0.40831610082646835, -0.00016188993599345562, -2.0011876919882003, -0.0764803841657652, 2.430139631960127, 0.8534804977636645] #~26cm distance
   #joint_ori = [-0.06033718608193325, 0.1957925676774354, 0.1446464792309258, -2.1242161722067974, -0.12300981136857973, 2.3759525292393855, -0.6157846782301644] #rotated ee grip
   franka = Franka(topic='/franka/', node_name='franka2_3_talker')
   #franka = Franka(init_node=True)
   franka.rate.sleep()

   #Import traj and duration from CSV
   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')

   traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_demoDMP.csv", delimiter=',') #NOTE: set name here!
   vel_traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_vel_demoDMP.csv", delimiter=',') #NOTE: set name here!

   dt = 0.001 #1/30 #/ 120 #NOTE: modify this to match generated CSV or set this FPS when generating csv!
   #NOTE: seems like I have to use higher dt than 1/FPS, even when there is a sleep in franka.py
   #NOTE: 10 FPS seems too few points so it is too jumpy!

   input("Move robots to origin")

   joint_ori = traj[0]
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin


   franka.close_grippers_middle()
   input("Close grippers")
   franka.close_grippers() #NEW


   tf = traj.shape[0] * dt

   print('tf:', tf)



   #TODO: make sure xyz axes of robots match with xyz axes from demo!

   
   #TODO: remove sleep times etc

   filepath = os.path.join(datafolder+"/"+'executed_joint_trajectory.csv')

   if os.path.exists(filepath):
      os.remove(filepath)

   input("Perform dynamic primitive")
   time.sleep(10) #sleep 10s when operating robots alone

   franka.move(move_type='jvt',params=vel_traj, traj_duration=tf)

   time.sleep(tf) #let motion finish before plotting and closing grippers

   real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!

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

   franka.open_grippers_middle()
   #franka.open_grippers()
   #franka.release_grippers() #NEW
   #franka.move(move_type='d',params=traj, traj_duration=4.0)
