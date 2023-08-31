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

   joint_ori = [0.00987823605505126, 0.374127258697815, 0.0510372705978866, -1.88723720136350, -0.326075606836512, 2.22057041702639, 0.261161445783366] #First setpoint in BagFlip csv (without -1), ORIGINAL GRIP DIRECTION

   #joint_ori = [-0.151392100200986, 0.374799255536351, 0.103110341127385, -1.89092065937875, 0.256956160347786, 2.23422794073271, -0.216915860033620] #First setpoint in BagFlip csv (WITH -1), ORIGINAL GRIP DIRECTION
   

   #joint_ori = [-0.4086566598369908, 0.06066233018972612, 0.4103449780140034, -2.305708428782938, -0.14869967789150004, 2.406735659296419, -0.6797833158672728] #y aligned with base x
   #joint_ori = [-0.31192760009611187, 0.018771605972182755, 0.4549330232667995, -2.3624764677273156, 0.06496027958723219, 2.431859529764697, -0.7104102848795065]
   
   #franka = Franka(topic='/franka/', node_name='franka2_3_talker')
   franka = Franka(init_node=True)
   franka.rate.sleep()

   input("Move robots to origin")
   
   #pose = np.array([0.57474, -0.04, 0.20, 0.99882, -0.039110, 0.010790, 0.026625]).reshape(1,7)
   #pose = np.array([0.60, 0, 0.20, 1, 0, 0, 0]).reshape(1,7)
   #pose = np.array([0.75, 0, 0.20, 0.70752, 0.70544, 0.00497, 0.04151]).reshape(1,7) #new ee orientation
   #print("pose x:", pose[0])
   
   #first apply joint movement if it is far from desired location (large motion not doable in single linear motion)
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0) #for joint movement to origin
   franka.close_grippers() #NEW
   
   #then apply linear motion to get to correct distance to other robot depending on bag weight
   #franka.move(move_type='o', params=pose, traj_duration=1) #for linear movement to origin - didn't work with single point? Use DMP code to set offset?

   # Defined as X, Y, Z
   #euler_init = [-177.31653135,   -4.27029547,   -4.35228107]  # Before = [-179.24297822, -2.93434324, 42.65659799]


   dt = 0.001 #1/30 #/ 120 #NOTE: modify this to match generated CSV or set this FPS when generating csv!
   #NOTE: seems like I have to use higher dt than 1/FPS, even when there is a sleep in franka.py


   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')

   #Import traj and duration from CSV


   traj = np.genfromtxt(datafolder+"/trajectories/"+"joint_BagFlip.csv", delimiter=',') #NOTE: set name here!
   #NOTE: 10 FPS seems too few points so it is too jumpy!

   #traj = traj[0,:].reshape(1, 7)
   #traj[:,2] += 0.15
   #print("traj shape:", traj.shape)
   #NOTE: robot handles quaternions (X,Y,Z,W) https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.as_quat.html 
   # and traj should be (x,y,z,) >>> DMPfunc.py now generates in this order

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

   input("Perform dynamic primitive")

   #for i in range(2500, traj.shape[0]):
   #   #print("shape:", traj[i,:].reshape(1,-1).shape)
   #   print("i:", i)
   #   franka.movel(traj[i,:].reshape(1,-1), traj_duration=dt)

   # def quaternion_multiply(quaternion1, quaternion0):
   #    x0, y0, z0, w0 = quaternion0
   #    x1, y1, z1, w1 = quaternion1
   #    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
   #                      x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
   #                      -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
   #                      x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)


   #q_rot_y = [0, 1, 0, 0]
   #q_rot_newx = [-0.3420201, 0, 0, 0.9396926]

   #traj = traj[0,:].reshape(1, 7)
   
   #traj[:,3] *= -1
   #traj[:,6] *= -1

   #traj = np.array([[0.57548, -0.0978, 0.24847, 0.878821, -0.0368, -0.02587, 0.475]]) #easy case and latest with origin

   #traj = np.array([[0.559, 0.088, 0.2536, 0.901, -0.0089, -0.0012, -0.4331]]) #hard case and latest with origin


   #traj = np.array([[0.572, 0.023, 0.209, 0.6749, 0.6607, -0.23487, 0.22982]]) #tries to rotate all the way around

   # r = R.from_quat([0.901, -0.0089, -0.0012, -0.4331])
   # angles = r.as_euler('xyz', degrees=False)
   # rot_x = angles[0]
   # print("rot x:", rot_x)
   # rot_y = angles[1]
   # rot_z = angles[2]

   # print("rot y:", rot_y)
   # if rot_x < 0:
   #       # if rot_y > 0:
   #       #    rot_x = -3.14
   #       #    quat = (R.from_euler('XYZ', [rot_x,rot_y,0], degrees=False).as_quat())
   #       # else:
   #          rot_y = -3.14
   #          #rot_x = -3.14-rot_x
   #          rot_x = rot_x+3.14

   #          quat = (R.from_euler('XYZ', [rot_x,rot_y,rot_z], degrees=False).as_quat())
   # else:
   #    quat = (R.from_euler('XYZ', [rot_x,rot_y,rot_z], degrees=False).as_quat())

   
   # traj[:,3], traj[:,4], traj[:,5], traj[:,6] = quat[0], quat[1], quat[2], quat[3]

   # print(rot_x, rot_y, " -- quat: ",quat)


   #traj = np.array([[0.572, 0.023, 0.209, 0.6547, 0.63439, 0.28806, -0.29317]]) #

   #traj[:,2] += 0.15
   #traj[:,3] =  -0.3420201
   #traj[:,4] =  0.940
   #traj[:,6] =   0.342

   #q_orig = [traj[:,3], traj[:,4], traj[:,5], traj[:,6]]
   #q_new1 = quaternion_multiply(q_orig, q_rot_y)
   #q_new2 = quaternion_multiply(q_new1, q_rot_newx)
   #traj[:,3], traj[:,4], traj[:,5], traj[:,6] = q_new2[0], q_new2[1], q_new2[2], q_new2[3]
   #print("q_new2:", q_new2)
   
   #franka.movel(traj, traj_duration=dt)
   franka.move(move_type='jt',params=traj, traj_duration=tf)

   time.sleep(tf) #let motion finish before plotting and closing grippers

   real_traj = np.genfromtxt(filepath, delimiter=',') #NOTE: set name here!

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

   franka.open_grippers()
   #franka.release_grippers() #NEW
   #franka.move(move_type='d',params=traj, traj_duration=4.0)
