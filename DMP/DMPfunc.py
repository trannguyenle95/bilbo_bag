"""
Copyright (C) 2016 Travis DeWolf

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

#If using unedited
#import pydmps
#print(pydmps.__file__) #If I want to see where actual module is stored
#import pydmps.dmp_discrete

import ownDMP

def DMPfunc(filename, main_axis, inFPS = 120, outFPS=120, plot=False, save=True):
    #load OptriTrack csv and preprocess

    df = pd.read_csv("Demos/"+filename, skiprows = 6) #keep 7th row for df headers
    df = df.ffill(axis ='rows').bfill(axis ='rows') #copy values if there are gaps, so quaternions should still be valid (in Matlab I did interpolation)
    data = df.to_numpy()

    #assume start in origin plus offset
    data[:,6] = data[:,6] - data[0,6] 
    data[:,7] = data[:,7] - data[0,7]
    data[:,8] = data[:,8] - data[0,8]
    
    max_h = np.max(data[:,7])
    scale = 1
    if max_h > 0.6:
        scale = 0.6/max_h

    data[:,6:9] *= scale

    run_time = data[-1, 1]
    main_axis = 'x' #main rotation axis
    if main_axis == 'x':
        axis = 2
        data[:,2+4] = 0 #set linear movement in this axis to zero
    elif main_axis == 'y':
        axis = 3
        data[:,3+4] = 0 #set linear movement in this axis to zero
    elif main_axis == 'z':
        axis = 4
        data[:,4+4] = 0 #set linear movement in this axis to zero
    else: #assume x otherwise by default
        raise Exception("Unrecognized main axis to preprocess data with.")

    q_ax = data[:, axis].copy()
    qw_new = np.sqrt(1-q_ax**2) * np.sign(data[:, 5])
    
    theta = 2*np.arctan2(np.sqrt(q_ax**2), qw_new)
    q_ax = np.sin((theta - np.pi)/2) 
    qw_new = np.cos((theta - np.pi)/2) * (-1) #-1 for movement towards the camera with the robot...

    #TODO: keep old quaternion but rotate it using some function e.g. scipy
    #theta = 2*np.atan2(np.sqrt(data[:, axis]**2), data[:, 2])
    #q_ax = np.sin((theta + np.pi)/2)
    #qw_new = np.cos(np.arcsin(q_ax) + np.pi) * np.sign(data[:, 5] /2)
    #qw_new = np.cos((theta + np.pi)/2) # * np.sign(data[:, 5] /2)
    
    
    #Also change order of columns to q = [qw qx qy qz] here!
    data[:, 2:6] = 0
    data[:, axis+1] = q_ax
    data[:, 2] = qw_new

    dims = 7 #(qw, qx, qy, qz) and (x_pos, y_pos, z_pos)
    y_demo = data[:,2:2+dims].T

    #Init discrete DMP
    dmp = ownDMP.dmp_discrete.DMPs_discrete(n_dmps=dims, n_bfs=500, dt=1/inFPS, ay=np.ones(dims) * 10.0, run_time = run_time) #NOTE: modified cs so now runtime can be given as input in kwargs

    #Train from demo
    dmp.imitate_path(y_des=y_demo, plot=False)

    #Generate new trajectory from DMP
    y_track, dy_track, ddy_track = dmp.rollout(tau=inFPS/outFPS)

    #TODO: NORMALIZE QUATERNION OUTPUTS WHEN LEARNING EACH DIM INDIVIDUALLY! >> Should be done inside function before creating y_dt and y_ddt!
    #If it is enough for frankas to get position and angle (no velocities or acclerations) then normalization can maybe be done here?
    y_track[:,0:4] = y_track[:,0:4] / np.linalg.norm(y_track[:,0:4], axis = 1, keepdims = True)

    reorder = np.zeros_like(y_track) #change order to match (x,y,z,qx,qy,qz,qw) used by robot script
    reorder[:,0] = y_track[:, 4] #x
    reorder[:,1] = y_track[:, 6] * -1 #y = old z flipped
    reorder[:,2] = y_track[:, 5] #z = old y
    reorder[:,3] = y_track[:, 1] #qx
    reorder[:,4] = y_track[:, 3] * -1 #qy = old qz flipped
    reorder[:,5] = y_track[:, 2] #qz = old qy
    reorder[:,6] = y_track[:, 0] #qw

    #make qw always be positive
    #qw_sign = np.sign(reorder[:,6])
    #reorder[:,3:7] = reorder[:,3:7] * qw_sign.reshape(-1, 1)

    #add offsets
    reorder[:,0] += 0.575
    reorder[:,1]
    reorder[:,2] += 0.20

    #Plotting
    if plot == True:
        plt.figure(1)
        plt.plot(reorder[:, 0], label="pos_x")
        plt.plot(reorder[:, 1], label="pos_y")
        plt.plot(reorder[:, 2], label="pos_z")
        plt.title("Position components")
        plt.plot(data[:, 6], label="pos_x in")
        plt.plot(data[:, 7], label="pos_y in")
        plt.plot(data[:, 8], label="pos_z in")
        plt.legend()
    
        plt.figure(2)
        plt.plot(reorder[:, 3], label="qx")
        plt.plot(reorder[:, 4], label="qy")
        plt.plot(reorder[:, 5], label="qz")
        plt.plot(reorder[:, 6], label="qw")
        plt.title("Quaternion components")
        plt.plot(data[:, axis+1], label="q_ax in")
        plt.plot(data[:, 2], label="q_w in")
        plt.legend()

        plt.show()

    #save DMP to file
    path = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'DMP', 'DMP_'+filename)
    np.savetxt(path, reorder, delimiter=",")


if __name__ == '__main__':
    DMPfunc("sack_from_bag2.csv", "x", inFPS = 120, outFPS=1000, plot=True, save=True) #1000 Hz playback times factor to slow it down by