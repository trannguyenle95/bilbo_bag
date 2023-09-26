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
   #^above used when initial ee grip is not rotated
   #joint_ori = [-0.06033718608193325, 0.1957925676774354, 0.1446464792309258, -2.1242161722067974, -0.12300981136857973, 2.3759525292393855, -0.6157846782301644] #rotated ee grip

   #joint_ori = [-0.001637353693321636, -0.7838729663498644, 0.001332866853098354, -2.354591992496418, -0.004458093395537303, 1.5761728843787417, 0.775939115057535] #setpoint in move_to_start Franka ROS
   
   #joint_ori = [-0.0689248124066240, 0.242700503048464, 0.131427626426893, -2.08512240051222, -0.359692141974038, 2.27848131322344, -1.27023685738111] #First setpoint in BagFlip csv (wihtout -1)
   #joint_ori = [-0.447564825600187, 0.258926980573075, 0.403231110458737, -2.09391016585431, 0.188238139907825, 2.31048277296264, -1.75515831222447] #First setpoint in BagFlip csv (WITH -1)

   #joint_ori = [0.00987823605505126, 0.374127258697815, 0.0510372705978866, -1.88723720136350, -0.326075606836512, 2.22057041702639, 0.261161445783366] #First setpoint in BagFlip csv (without -1), ORIGINAL GRIP DIRECTION

   #joint_ori = [-0.151392100200986, 0.374799255536351, 0.103110341127385, -1.89092065937875, 0.256956160347786, 2.23422794073271, -0.216915860033620] #First setpoint in BagFlip csv (WITH -1), ORIGINAL GRIP DIRECTION
   
   #joint_ori = [0, 0.2837448589662732, 0, -2.0720574669683027, 0, 2.405712411822974, 0.7542077567525343] #NOTE: position [0.59,0,0.20] gripper in original orientation

   #joint_ori = [0.0711714439885711, 0.373619455135678, -0.0160332506422611, -1.88868112924497, -0.296405612466853, 2.22703550315137, 0.241982022300288]
   
   #joint_ori = [-0.001637353693321636, -0.7838729663498644, 0.001332866853098354, -2.354591992496418, -0.004458093395537303, 1.5761728843787417, 0.775939115057535] #Init position before delta_t updates
   
   #joint_ori = [-0.4086566598369908, 0.06066233018972612, 0.4103449780140034, -2.305708428782938, -0.14869967789150004, 2.406735659296419, -0.6797833158672728] #y aligned with base x
   #joint_ori = [-0.31192760009611187, 0.018771605972182755, 0.4549330232667995, -2.3624764677273156, 0.06496027958723219, 2.431859529764697, -0.7104102848795065]
   
   #franka = Franka(topic='/franka/', node_name='franka2_3_talker')
   franka = Franka(init_node=True)
   franka.rate.sleep()

   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')

   #Import traj and duration from CSV
   traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_2h_flip2.csv", delimiter=',') #NOTE: set name here!
   dt = 0.001 #1/30 #/ 120 #NOTE: modify this to match generated CSV or set this FPS when generating csv!
   #NOTE: seems like I have to use higher dt than 1/FPS, even when there is a sleep in franka.py

   input("Move robots to origin")
   
   joint_ori = traj[0]

   #first apply joint movement if it is far from desired location (large motion not doable in single linear motion)
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin


   franka.close_grippers_middle()
   input("Close grippers")
   franka.close_grippers() #NEW
   
   '''
   #INTERPOLATION
   #https://stackoverflow.com/questions/16388110/double-the-length-of-a-python-numpy-array-with-interpolated-new-values
   print("traj shape:", traj.shape)

   N = len(traj[:,0])
   X = np.arange(0, 10*N, 10)
   X_new = np.arange(10*N-1)       # Where you want to interpolate
   j1_new = np.interp(X_new, X, traj[:,0]) 
   j2_new = np.interp(X_new, X, traj[:,1])
   j3_new = np.interp(X_new, X, traj[:,2]) 
   j4_new = np.interp(X_new, X, traj[:,3]) 
   j5_new = np.interp(X_new, X, traj[:,4]) 
   j6_new = np.interp(X_new, X, traj[:,5]) 
   j7_new = np.interp(X_new, X, traj[:,6]) 


   new_traj = np.vstack((j1_new, j2_new, j3_new, j4_new, j5_new, j6_new, j7_new)).T

   print("new traj shape:", new_traj.shape)

   plt.figure(2)
   plt.plot(traj[:, 0], 'r--', label="j1")
   plt.plot(traj[:, 1],'g--', label="j2")
   plt.plot(traj[:, 2], 'b--', label="j3")
   plt.show()
   plt.figure(3)
   plt.plot(new_traj[:, 0], 'r-', label="new_j1")
   plt.plot(new_traj[:, 1],'g-', label="new_j2")
   plt.plot(new_traj[:, 2], 'b-', label="new_j3")
   plt.legend()
   plt.show()

   plt.figure(4)
   plt.plot(traj[:, 3], 'r--', label="j4")
   plt.plot(traj[:, 4],'g--', label="j5")
   plt.plot(traj[:, 5], 'b--', label="j6")
   plt.plot(traj[:, 6], 'm--', label="j7")
   plt.show()
   plt.figure(5)
   plt.plot(new_traj[:, 3], 'r-', label="new_j1")
   plt.plot(new_traj[:, 4],'g-', label="new_j2")
   plt.plot(new_traj[:, 5], 'b-', label="new_j3")
   plt.plot(new_traj[:, 6], 'm-', label="new_j7")
   plt.legend()
   plt.show()

   traj = new_traj #NOTE: use interpolation
   '''
   #NOTE: using initial position as target each time
   #traj = np.repeat(np.array([[-0.06033718608193325, 0.1957925676774354, 0.1446464792309258, -2.1242161722067974, -0.12300981136857973, 2.3759525292393855, -0.6157846782301644]]), traj.shape[0], axis=0)
   #print("new traj shape:", traj.shape)

   #NOTE: using similar target as in Franka ROS example joint controller
   # t = np.arange(0,10,1/1000)
   # delta_angle = np.double(np.pi / 8 * (1 - np.cos(np.pi / 2.5 * t)))

   # print("t shape:", t.shape[0])
   # #print("shape delta:", delta_angle.shape)

   # traj = np.repeat(np.array([[-0.001637353693321636, -0.7838729663498644, 0.001332866853098354, -2.354591992496418, -0.004458093395537303, 1.5761728843787417, 0.775939115057535]], dtype="double"), t.shape[0], axis=0)

   # #traj[:,0] += delta_angle
   # #traj[:,1] += delta_angle
   # #traj[:,2] += delta_angle
   # traj[:,3] += delta_angle * 1
   # traj[:,4] += delta_angle * 1
   # #traj[:,5] += delta_angle
   # traj[:,6] += delta_angle * 1

   # plt.figure(4)
   # plt.plot(traj[:, 0], label="j1")
   # plt.plot(traj[:, 1], label="j2")
   # plt.plot(traj[:, 2], label="j3")
   # plt.plot(traj[:, 3], label="j4")
   # plt.plot(traj[:, 4], 'r--', label="j5")
   # plt.plot(traj[:, 5], label="j6")
   # plt.plot(traj[:, 6], label="j7")
   # plt.legend()
   # plt.show()

   # # traj[:,0] += delta_angle
   # # traj[:,1] += delta_angle
   # # traj[:,2] += delta_angle
   # # traj[:,3] += delta_angle
   # # traj[:,5] += delta_angle   
   # # traj[:,6] += delta_angle

   # print("shape traj:", traj.shape)

   tf = traj.shape[0] * dt

   print('tf:', tf)

   # plt.figure(1)
   # plt.plot(traj[:, 0], label="pos_x")
   # plt.plot(traj[:, 1], label="pos_y")
   # plt.plot(traj[:, 2], label="pos_z")
   # plt.title("Position components")
   # plt.legend()

   # plt.figure(2)
   # plt.plot(traj[:, 3], label="qx")
   # plt.plot(traj[:, 4], label="qy")
   # plt.plot(traj[:, 5], label="qz")
   # plt.plot(traj[:, 6], label="qw")
   # plt.title("Quaternion components")
   # plt.legend()

   # plt.show()


   #TODO: make sure xyz axes of robots match with xyz axes from demo!

   

   #TODO: remove sleep times etc

   #franka.movedynamic_ori(quintic_traj=quintic_traj, tf=tf)
   #franka.movedynamic_ori(quintic_traj=traj, tf=tf)

   filepath = os.path.join(datafolder+"/"+"executed_joint_trajectory.csv") #NOTE: new file for joint traj!

   if os.path.exists(filepath):
      os.remove(filepath) 


   pose_file = os.path.join(datafolder+"/"+"joint_control_pose.csv")
   if os.path.exists(pose_file):
      os.remove(pose_file) 


   joint_vel_file = os.path.join(datafolder+"/"+"executed_joint_velocities.csv")
   if os.path.exists(joint_vel_file):
      os.remove(joint_vel_file) 


   input("Perform dynamic primitive")


   franka.move(move_type='jt',params=traj, traj_duration=tf)

   time.sleep(tf) #let motion finish before plotting and closing grippers

   real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!


   #Plot pose from executed actual motion
   ref_pose = np.genfromtxt(datafolder+"/trajectories/"+"pose_2h_flip2.csv", delimiter=',') #NOTE: set name here!
   real_pose = np.genfromtxt(pose_file, delimiter=',') #NOTE: set name here!
   real_vel = np.genfromtxt(joint_vel_file, delimiter=',') #NOTE: set name here!
   
   real_acc = np.diff(real_vel) / 0.001; 
   print("real acc shape:", real_acc.shape)

   #TODO: add more plots
   plt.figure(1)
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

   #PLOT Cartesian Pose during joint control 
   plt.figure(2)
   plt.plot(ref_pose[:, 0], 'r-', label="ref pos_x")
   plt.plot(ref_pose[:, 1],'g-', label="ref pos_y")
   plt.plot(ref_pose[:, 2], 'b-', label="ref pos_z")
   plt.title("Position components")
   plt.plot(real_pose[:, 0], 'r--', label="actual pos_x")
   plt.plot(real_pose[:, 1], 'g--', label="actual pos_y")
   plt.plot(real_pose[:, 2], 'b--', label="actual pos_z")
   plt.legend()

   plt.figure(3)
   plt.plot(ref_pose[:, 3], 'r-', label="ref qx")
   plt.plot(ref_pose[:, 4], 'g-', label="ref qy")
   plt.plot(ref_pose[:, 5], 'b-', label="ref qz")
   plt.plot(ref_pose[:, 6], 'm-', label="ref qw")
   plt.title("Quaternion components")
   plt.plot(real_pose[:, 3], 'r--', label="actual qx")
   plt.plot(real_pose[:, 4], 'g--', label="actual qy")
   plt.plot(real_pose[:, 5], 'b--', label="actual qz")
   plt.plot(real_pose[:, 6], 'm--', label="actual qw")
   plt.legend()

   plt.figure(4)
   plt.title("Joint vels")
   plt.plot(real_vel[:, 0], '-', label="vel j1")
   plt.plot(real_vel[:, 1],'-', label="vel j2")
   plt.plot(real_vel[:, 2], '-', label="vel j3")
   plt.plot(real_vel[:, 3], '-', label="vel j4")
   plt.plot(real_vel[:, 4], '-', label="vel j5")
   plt.plot(real_vel[:, 5], '-', label="vel j6")
   plt.plot(real_vel[:, 6], '-', label="vel j7")
   plt.legend()
   plt.show()

   plt.figure(5)
   plt.title("Joint accs")
   plt.plot(real_acc[:, 0], '-', label="acc j1")
   plt.plot(real_acc[:, 1],'-', label="acc j2")
   plt.plot(real_acc[:, 2], '-', label="acc j3")
   plt.plot(real_acc[:, 3], '-', label="acc j4")
   plt.plot(real_acc[:, 4], '-', label="acc j5")
   plt.plot(real_acc[:, 5], '-', label="acc j6")
   plt.plot(real_acc[:, 6], '-', label="acc j7")
   plt.legend()
   plt.show()



   plt.show()

   franka.open_grippers()
   #franka.release_grippers() #NEW
   #franka.move(move_type='d',params=traj, traj_duration=4.0)
