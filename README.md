# Dynamic Manipulation of Deformable Objects using Imitation Learning with Adaptation to Hardware Constraints
This page contains code from the following project: https://sites.google.com/view/bilbo-bag.
Below is a brief description of the contents of each folder:

**Data** contais subfolders for human demonstrations, the recorded bag states during experiments, and robot trajectories of the most recent run (mainly used for generating plots for debugging purposes).

**SupportScripts** contains MATLAB and Python scripts for processing human demonstrations so that they can be used in the constrained DMPs, for estimating the state of the bag from a recorded cloud of marker positions, and for generating plots.

**TC_DMP_constrainedVelAcc** contains code adapted from https://github.com/albindgit/TC_DMP_constrainedVelAcc, which is used to implement two constrained DMP versions.

**franka_david** contains code for running the robots. The subfolder "scripts" contains the logic of the experiments, with the exception of "franka.py". This file and the remaining folders contain code related to implementing the controllers and motion-generators for the robots, and it builds upon code used in previous projects of the research group. As the robot control code is not the key focus of this project we do not prioritize its readability. 

**novel-DMP-constraints-main** contains code adapted from https://github.com/Slifer64/novel-DMP-constraints, which is adapted to implement a third constrained DMP version.

**skymul_natnet_ros_cpp**: we use the following ROS driver used to record the markers on the bag https://github.com/SkyMul/skymul_natnet_ros_cpp/tree/main.
