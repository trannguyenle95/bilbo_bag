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

import pydmps
import pydmps.dmp_discrete


#load OptriTrack csv and preprocess

df = pd.read_csv("BagFlip.csv", skiprows = 6) #keep 7th row for df headers
df = df.ffill(axis ='rows').bfill(axis ='rows') #copy values if there are gaps, so quaternions should still be valid (in Matlab I did interpolation)
data = df.to_numpy()


demo_timesteps = data.shape[0]
#print("demo t:", demo_timesteps)

main_axis = 'x' #main rotation axis

if main_axis == 'x':
    axis = 2
elif main_axis == 'y':
    axis = 3
elif main_axis == 'z':
    axis = 4
else: #assume x otherwise by default
    raise Exception("Unrecognized main axis to preprocess data with.")

q_ax = data[:, axis].copy()
qw_new = np.sqrt(1-q_ax**2) * np.sign(data[:, 5])
#Also change order of columns to q = [qw qx qy qz] here!
data[:, 2:6] = 0
data[:, axis+1] = q_ax
data[:, 2] = qw_new


dims = 7 #(qw, qx, qy, qz) and (x_pos, y_pos, z_pos)
y_demo = data[:,2:2+dims].T

# test normal run
dmp = pydmps.dmp_discrete.DMPs_discrete(n_dmps=dims, n_bfs=500, ay=np.ones(dims) * 10.0)
#y_track = []
#dy_track = []
#ddy_track = []

dmp.imitate_path(y_des=y_demo, plot=False)
y_track, dy_track, ddy_track = dmp.rollout(tau=100/demo_timesteps) #default timesteps are 100 with tau=1

#print("output shape:", y_track.shape)

#TODO: NORMALIZE QUATERNION OUTPUTS WHEN LEARNING EACH DIM INDIVIDUALLY! >> Should be done inside function before creating y_dt and y_ddt!
#If it is enough for frankas to get position and angle (no velocities or acclerations) then normalization can maybe be done here?
y_track[:,0:4] = y_track[:,0:4] / np.linalg.norm(y_track[:,0:4], axis = 1, keepdims = True)



#TODO: how does output shape from discrete dmp work > what kind of scaling is applied and how to alter it??
# DMP sets runtime to 1s, which can then be scaled with tau in output... (defualt 100 step output trajectory)


#Plotting
plt.figure(1)
plt.plot(y_track[:, 0], label="qw")
plt.plot(y_track[:, 1], label="qx")
plt.plot(y_track[:, 2], label="qy")
plt.plot(y_track[:, 3], label="qz")
plt.title("Quaternion components")
plt.plot(data[:, axis+1], label="q_ax in")
plt.plot(data[:, 2], label="q_w in")
plt.legend()

plt.figure(2)
plt.plot(y_track[:, 4], label="pos_x")
plt.plot(y_track[:, 5], label="pos_y")
plt.plot(y_track[:, 6], label="pos_z")
plt.title("Position components")
plt.plot(data[:, 6], label="pos_x in")
plt.plot(data[:, 7], label="pos_y in")
plt.plot(data[:, 8], label="pos_z in")

plt.legend()

plt.show()


#save DMP to file
#np.savetxt("outputDMP.csv", y_track, delimiter=",")