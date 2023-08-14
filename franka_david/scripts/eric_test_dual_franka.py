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
   joint_ori = [-0.04978088093284428, 0.40831610082646835, -0.00016188993599345562, -2.0011876919882003, -0.0764803841657652, 2.430139631960127, 0.8534804977636645] #~26cm distance
   franka = Franka(topic='/franka/', node_name='franka2_3_talker')
   #franka = Franka(init_node=True)
   franka.rate.sleep()

   input("Move robots to origin")
   franka.move(move_type='j', params=joint_ori, traj_duration=3.0)

   # Defined as X, Y, Z
   #euler_init = [-177.31653135,   -4.27029547,   -4.35228107]  # Before = [-179.24297822, -2.93434324, 42.65659799]


   dt = 0.001 #1/30 #/ 120 #NOTE: modify this to match generated CSV or set this FPS when generating csv!
   #NOTE: seems like I have to use higher dt than 1/FPS, even when there is a sleep in franka.py


   #Import traj and duration from CSV
   datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'DMP')

   traj = np.genfromtxt(datafolder+"/DMP/"+"DMP_sack_from_bag2.csv", delimiter=',') #NOTE: set name here!
   #NOTE: 10 FPS seems too few points so it is too jumpy!

   #traj = traj[0,:].reshape(1, 7)
   #traj[:,2] += 0.15
   #print("traj shape:", traj.shape)
   #NOTE: robot handles quaternions (X,Y,Z,W) https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.as_quat.html 
   # and traj should be (x,y,z,) >>> DMPfunc.py now generates in this order
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
   filepath = os.path.join(datafolder+"/"+'executed_trajectory.csv')

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

   plt.show()

    #franka.move(move_type='d',params=traj, traj_duration=4.0)
