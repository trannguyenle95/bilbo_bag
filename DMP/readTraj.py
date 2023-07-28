import numpy as np
import matplotlib.pyplot as plt

def readTraj(filename, fps):
    traj = np.genfromtxt(filename, delimiter=',') #NOTE: set name here!
    tf = traj.shape[0] / fps

    return traj, tf

if __name__ == '__main__':
    traj, tf = readTraj("120fps_DMP_BagFlip.csv", fps = 1000)
    print('traj shape:', traj.shape)
    print('dt:', tf)

    plt.figure(1)
    plt.plot(traj[:, 0], label="pos_x")
    plt.plot(traj[:, 1], label="pos_y")
    plt.plot(traj[:, 2], label="pos_z")
    plt.title("Position components")
    plt.legend()

    plt.figure(2)
    plt.plot(traj[:, 3], label="qx")
    plt.plot(traj[:, 4], label="qy")
    plt.plot(traj[:, 5], label="qz")
    plt.plot(traj[:, 6], label="qw")
    plt.title("Quaternion components")
    plt.legend()

    plt.show()
