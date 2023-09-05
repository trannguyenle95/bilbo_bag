#pragma once

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <controller_interface/controller_base.h>
#include <geometry_msgs/PoseStamped.h>
#include <ros/ros.h>


#include <franka/robot.h>
#include <franka/gripper.h>

#include <matrix.hpp>
#include <Frame3D.h>

#include <franka_david/MotionPy.h>
#include <franka_david/MovetoPy.h>
#include <franka_david/GripperPy.h>
#include <franka_david/GripperGraspPy.h>
#include <franka_david/JointMotionPy.h>
#include <franka_david/JointTrajPy.h>

class CartesianPythonController
{
public:

    // Subscriber for motion
    ros::Subscriber _sub_motion;
    ros::Subscriber _sub_motion_ori;
    // Subscriber for moveto
    ros::Subscriber _sub_moveto;
    // Subscriber joint motion
    ros::Subscriber _sub_joint;
    // Subscriber joint trajectory
    ros::Subscriber _sub_joint_traj;
    // Subscriber gripper
    ros::Subscriber _sub_gripper_move;
    // Subscriber gripper
    ros::Subscriber _sub_gripper_grasp;
    
    // Publisher for reached state
    ros::Publisher _pub_reached;
    
    // Goal pose
    Eigen::Vector3d _goal_pose;
    // Current pose
    Eigen::Vector3d _current_pose;
    
    // Delta time - 1KHz control
    double _dt = 0.001;
    
    // Method that subscribes to a joint motion
    void jointMotionCallback(const franka_david::JointMotionPyPtr& msg);
    
    // Method that subscribes to a joint trajectory, ADDED 29.08
    void jointTrajectoryCallback(const franka_david::JointTrajPyPtr& msg);
    
    // Method that subscribes to motion
    void motionCallback(const franka_david::MotionPyPtr& msg);
    void motionOriCallback(const franka_david::MotionPyPtr& msg);
    // Method that subscribes to moveto
    void movetoCallback(const franka_david::MovetoPyPtr& msg);
    
    // Method that subscribes to gripper
    void gripperMoveCallback(const franka_david::GripperPyPtr& msg);
    
    // Method that subscribes to gripper
    void gripperGraspCallback(const franka_david::GripperGraspPyPtr& msg);
    
    
    CartesianPythonController();
    
    void runControl(math::Transform3D* trajectory, int N);
        
// private:
    franka::Robot* _robot;
    franka::Gripper* _gripper;
    bool _franka3 = false; //change to either 2 or 3 to match the robot
};
