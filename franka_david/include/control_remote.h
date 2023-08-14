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
#include <franka_david/RemoteMotionPy.h>


class CartesianRemoteController
{
public:

    // Subscriber for motion
    ros::Subscriber _sub_run_motion;
    
     // Subscriber for motion
    ros::Subscriber _sub_moveto;
    
    // Subscriber for trajetory motion
    ros::Subscriber _sub_motion;
    
    // Subscriber for trajectory motion with orientation
    ros::Subscriber _sub_motion_ori;
    
    // Subscriber for joint
    ros::Subscriber _sub_joint;
    
    // Publisher for reached state
    ros::Publisher _pub_reached;
    
    // Goal pose
    Eigen::Vector3d _goal_pose;
    // Current pose
    Eigen::Vector3d _current_pose;
    
    // Delta time - 1KHz control
    double _dt = 0.001;
    
    // Method that subscribes to a joint motion
    void runRemoteMotionCallback(const franka_david::RemoteMotionPyPtr& msg);
    
    // Method that subscribes to moveto
    void movetoCallback(const franka_david::MovetoPyPtr& msg);
    
    // Method to execute full trajectory
    void motionCallback(const franka_david::MotionPyPtr& msg);
    
    // Method to execute full trajectory changing the orientation
    void motionOriCallback(const franka_david::MotionPyPtr& msg);
    
    // Method to execute joint motion
    void jointMotionCallback(const franka_david::JointMotionPyPtr& msg);
    
    
    CartesianRemoteController();
    
    void runControl(math::Transform3D* trajectory, int N);
        
// private:
    franka::Robot* _robot;
    franka::Gripper* _gripper;
    bool _franka3 = false; //change this to match current robot (2 or 3)
    
};
