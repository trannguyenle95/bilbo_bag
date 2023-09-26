#!/bin/sh
roslaunch franka_control franka_control.launch robot_ip:=172.16.0.2 & sleep 3;
rostopic pub -1 /franka_control/error_recovery/goal franka_msgs/ErrorRecoveryActionGoal "{}"
