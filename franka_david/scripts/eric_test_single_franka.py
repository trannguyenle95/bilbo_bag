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
   #rospy.init_node('franka3_talker', anonymous=True)
   
   #For getting joint states:
   #roslaunch franka_visualization franka_visualization.launch robot_ip:=172.16.0.2 load_gripper:=True
   #Then run:
   #rostopic echo joint_states
   
   #joint_ori = [-0.06153707569372121, 0.23268072162435294, -0.003733379824253959, -2.120620626949313, -0.07440938119840552, 2.374850448676014, 0.851590066155449] #hori x align
   #joint_ori = [-0.04978088093284428, 0.40831610082646835, -0.00016188993599345562, -2.0011876919882003, -0.0764803841657652, 2.430139631960127, 0.792066740804676] #~26cm distance

   joint_ori = [0, 0.2837448589662732, 0, -2.0720574669683027, 0, 2.405712411822974, 0.7542077567525343] #NOTE: position [0.59,0,0.20] gripper in original orientation

   #joint_ori = [0.0772848981560748, 0.240594409923399, -0.0202525725346743, -2.08784449998159, -0.31363716946651, 2.2895568232216, 0.269719648304103]

   #^above used when initial ee grip is not rotated
   #joint_ori = [-0.06033718608193325, 0.1957925676774354, 0.1446464792309258, -2.1242161722067974, -0.12300981136857973, 2.3759525292393855, -0.6157846782301644] #rotated ee grip
   
   #joint_ori = [-0.4086566598369908, 0.06066233018972612, 0.4103449780140034, -2.305708428782938, -0.14869967789150004, 2.406735659296419, -0.6797833158672728] #y aligned with base x
   #joint_ori = [-0.31192760009611187, 0.018771605972182755, 0.4549330232667995, -2.3624764677273156, 0.06496027958723219, 2.431859529764697, -0.7104102848795065]
   
   #franka = Franka(topic='/franka/', node_name='franka2_3_talker')


   franka = Franka(init_node=True)
   franka.rate.sleep()

   #Import traj and duration from CSV
   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')

   traj = np.genfromtxt(datafolder+"/trajectories/"+"pose_sack_from_bag1.csv", delimiter=',') #NOTE: set name here!
   #NOTE: 10 FPS seems too few points so it is too jumpy!

   
   input("Move robots to origin")
   
   dt = 0.001 #1/30 #/ 120 #NOTE: modify this to match generated CSV or set this FPS when generating csv!
   #NOTE: seems like I have to use higher dt than 1/FPS, even when there is a sleep in franka.py

   tf = traj.shape[0] * dt
   print('tf:', tf)

   #first apply joint movement if it is far from desired location (large motion not doable in single linear motion)
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
   franka.close_grippers() #NEW
   





   #TODO: make sure xyz axes of robots match with xyz axes from demo!

   #TODO: remove sleep times etc


   filepath = os.path.join(datafolder+"/"+"executed_trajectory.csv")
   if os.path.exists(filepath):
      os.remove(filepath) 

   joint_file = os.path.join(datafolder+"/"+"pose_control_joints.csv")
   if os.path.exists(joint_file):
      os.remove(joint_file) 

   input("Perform dynamic primitive")

   franka.move(move_type='o',params=traj, traj_duration=tf)

   real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!


   plt.figure(1)
   plt.plot(traj[:, 0], 'r-', label="pos_x")
   plt.plot(traj[:, 1],'g-', label="pos_y")
   plt.plot(traj[:, 2], 'b-', label="pos_z")
   plt.title("Position components")
   plt.plot(real_traj[:, 0], 'r--', label="pos_x")
   plt.plot(real_traj[:, 1], 'g--', label="pos_y")
   plt.plot(real_traj[:, 2], 'b--', label="pos_z")
   plt.legend()

   plt.figure(2)
   plt.plot(traj[:, 3], 'r-', label="qx")
   plt.plot(traj[:, 4], 'g-', label="qy")
   plt.plot(traj[:, 5], 'b-', label="qz")
   plt.plot(traj[:, 6], 'm-', label="qw")
   plt.title("Quaternion components")
   plt.plot(real_traj[:, 3], 'r--', label="qx")
   plt.plot(real_traj[:, 4], 'g--', label="qy")
   plt.plot(real_traj[:, 5], 'b--', label="qz")
   plt.plot(real_traj[:, 6], 'm--', label="qw")
   plt.legend()

   #Plot joints from executed actual motion
   ref_joints = np.genfromtxt(datafolder+"/trajectories/"+"joint_sack_from_bag1.csv", delimiter=',') #NOTE: set name here!
   real_joints = np.genfromtxt(joint_file, delimiter=',') #NOTE: set name here!

   plt.figure(3)
   plt.plot(ref_joints[:, 0], '-', label="ref j1")
   plt.plot(ref_joints[:, 1],'-', label="ref j2")
   plt.plot(ref_joints[:, 2], '-', label="refj3")
   plt.plot(ref_joints[:, 3], '-', label="ref j4")
   plt.plot(ref_joints[:, 4], '-', label="ref j5")
   plt.plot(ref_joints[:, 5], '-', label="ref j6")
   plt.plot(ref_joints[:, 6], '-', label="ref j7")
   plt.title("Real joint values")
   plt.plot(real_joints[:, 0], '--', label="j1")
   plt.plot(real_joints[:, 1], '--', label="j2")
   plt.plot(real_joints[:, 2], '--', label="j3")
   plt.plot(real_joints[:, 3], '--', label="j4")
   plt.plot(real_joints[:, 4], '--', label="j5")
   plt.plot(real_joints[:, 5], '--', label="j6")
   plt.plot(real_joints[:, 6], '--', label="j7")
   plt.legend()

   plt.figure(4)
   plt.title("joint 1")
   plt.plot(ref_joints[:, 0], '-', label="ref")
   plt.plot(real_joints[:, 0], '--', label="actual")
   plt.legend()

   plt.figure(5)
   plt.title("joint 2")
   plt.plot(ref_joints[:, 1], '-', label="ref")
   plt.plot(real_joints[:, 1], '--', label="actual")
   plt.legend()

   plt.figure(6)
   plt.title("joint 3")
   plt.plot(ref_joints[:, 2], '-', label="ref")
   plt.plot(real_joints[:, 2], '--', label="actual")
   plt.legend()

   plt.figure(7)
   plt.title("joint 4")
   plt.plot(ref_joints[:, 3], '-', label="ref")
   plt.plot(real_joints[:, 3], '--', label="actual")
   plt.legend()

   plt.figure(8)
   plt.title("joint 5")
   plt.plot(ref_joints[:, 4], '-', label="ref")
   plt.plot(real_joints[:, 4], '--', label="actual")
   plt.legend()

   plt.figure(9)
   plt.title("joint 6")
   plt.plot(ref_joints[:, 5], '-', label="ref")
   plt.plot(real_joints[:, 5], '--', label="actual")
   plt.legend()

   plt.figure(10)
   plt.title("joint 7")
   plt.plot(ref_joints[:, 6], '-', label="ref")
   plt.plot(real_joints[:, 6], '--', label="actual")
   plt.legend()


   plt.show()

   franka.open_grippers()
   #franka.release_grippers() #NEW
   #franka.move(move_type='d',params=traj, traj_duration=4.0)
