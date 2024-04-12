# dynamics_bag
This page contains code from the following project: https://sites.google.com/view/bilbo-bag.
Below is a brief description of the contents of each folder:

**Data** contais subfolders for human demonstrations, the recorded bag states during experiments, and robot trajectories of the most recent run (mainly used for generating plots for debugging purposes).
**SupportScripts** contains MATLAB and Python scripts for processing human demonstrations so that they can be used in the constrained DMPs, for estimating the state of the bag from a recorded cloud of marker positions, and for generating plots.
**TC_DMP_constrainedVelAcc** contains code adapted from https://github.com/albindgit/TC_DMP_constrainedVelAcc, which is used to implement two constrained DMP versions.
**franka_david** WIP
**novel-DMP-constraints-main** contains code adapted from https://github.com/Slifer64/novel-DMP-constraints, which is adapted to implement a third constrained DMP version.
**skymul_natnet_ros_cpp** WIP FIX LINK (links to the ROS driver used to record the markers for estimating the bag state)



#TODO:
- remove bagmetrics.sh as it is redundant now with the bagmetrics files in SupportScripts
- mention "running_the_robots" or move it into franka_david
- decide whether to include the franka_david folder for controlling robots?
