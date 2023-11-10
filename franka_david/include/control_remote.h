#pragma once

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>

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
#include <franka_david/JointTrajPy.h>

#include <franka_david/JointVelTrajPy.h>

#include <franka_david/MoveRelativePy.h>

#include <franka_david/DiffJointVelTrajPy.h>
#include <franka_david/DiffJointMotionPy.h>

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
    
    // Subscriber joint trajectory
    ros::Subscriber _sub_joint_traj;
    
    // Subscriber gripper
    ros::Subscriber _sub_gripper_move;
    // Subscriber gripper
    ros::Subscriber _sub_gripper_grasp;

    ros::Publisher _franka_2_state;
    ros::Publisher _franka_3_state;
    ros::Subscriber _sub_franka_2_state;
    ros::Subscriber _sub_franka_3_state;
    
    // Subscriber joint velocity trajectory added 26.09
    ros::Subscriber _sub_joint_vel_traj;
    
    // Subscriber with separate trajectories for each robot added 10.11
    ros::Subscriber _sub_diff_joint_vel_traj;
    ros::Subscriber _sub_diff_joint;
    
    // Subscriber relative position added 12.10
    ros::Subscriber _sub_move_relative;
    
    // Publisher for reached state
    ros::Publisher _pub_reached;
    
    // Goal pose
    Eigen::Vector3d _goal_pose;
    // Current pose
    Eigen::Vector3d _current_pose;
    
    // Delta time - 1KHz control
    double _dt = 0.001;

    // Method for publishing the state
    void publishFranka2State();
    // Method for publishing the state
    void publishFranka3State();
    
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
    
    // Method that subscribes to a joint trajectory, ADDED 05.09
    void jointTrajectoryCallback(const franka_david::JointTrajPyPtr& msg);
    
    // Method that subscribes to a joint velocity trajectory, ADDED 26.09
    void jointVelocityTrajectoryCallback(const franka_david::JointVelTrajPyPtr& msg);
    
    // separate joint velocity trajectory for franka2 and franka3
    void diffJointVelocityTrajectoryCallback(const franka_david::DiffJointVelTrajPyPtr& msg);
    void diffJointMotionCallback(const franka_david::DiffJointMotionPyPtr& msg);
    
    // Method that subscribes to gripper
    void gripperMoveCallback(const franka_david::GripperPyPtr& msg);
    
    // Method that subscribes to gripper
    void gripperGraspCallback(const franka_david::GripperGraspPyPtr& msg);
    
    // Method that subscribes to a relative motion, ADDED 12.10
    void moveRelativeCallback(const franka_david::MoveRelativePyPtr& msg);
    
    //Added error callback
    void franka2ErrorCallback(const std_msgs::BoolPtr& msg);
    void franka3ErrorCallback(const std_msgs::BoolPtr& msg);
    
    
    
    CartesianRemoteController();
    
    void runControl(math::Transform3D* trajectory, int N);
        
// private:
    franka::Robot* _robot;
    franka::Gripper* _gripper;
    bool _franka3 = false; //change this to match current robot (2 or 3)
        /// \brief Thread running the rosPublishQueue.
    std::thread _franka2statethread;
    std::thread _franka3statethread;
    
    
};
