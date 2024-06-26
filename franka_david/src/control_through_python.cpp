// Copyright (c) 2017 Franka Emika GmbH
// Use of this source code is governed by the Apache-2.0 license, see LICENSE

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <controller_interface/controller_base.h>
#include <franka_hw/franka_cartesian_command_interface.h>
#include <hardware_interface/hardware_interface.h>
#include <pluginlib/class_list_macros.h>

#include <geometry_msgs/PoseStamped.h>
#include <ros/ros.h>

#include <franka_hw/franka_cartesian_command_interface.h>
#include <franka/duration.h>
#include <franka/exception.h>
#include <franka/model.h>
#include <franka/robot.h>

#include <matrix.hpp>
#include <Frame3D.h>
#include <quaternion.hpp>
#include <examples_common.h>

#include <control_through_python.h>


// template <class T, size_t N>
// std::ostream& operator<<(std::ostream& ostream, const std::array<T, N>& array) {
//   ostream << "[";
//   std::copy(array.cbegin(), array.cend() - 1, std::ostream_iterator<T>(ostream, ","));
//   std::copy(array.cend() - 1, array.cend(), std::ostream_iterator<T>(ostream));
//   ostream << "]";
//   return ostream;
// }


CartesianPythonController::CartesianPythonController()
{   
    
    ros::NodeHandle control_py_node("py_control_node");
    this->_robot = new franka::Robot("172.16.0.2");
    this->_gripper = new franka::Gripper("172.16.0.2");
    
    
    // TODO Modify delay
    this->_sub_motion_ori = control_py_node.subscribe(
        "/frankapy/motion_ori", 20, &CartesianPythonController::motionOriCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());

    this->_sub_motion = control_py_node.subscribe(
        "/frankapy/motion", 20, &CartesianPythonController::motionCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());
    
    this->_sub_moveto = control_py_node.subscribe(
        "/frankapy/moveto", 20, &CartesianPythonController::movetoCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());
    
    this->_sub_gripper_move = control_py_node.subscribe(
        "/frankapy/gripper_move", 20, &CartesianPythonController::gripperMoveCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());
    
    this->_sub_gripper_grasp = control_py_node.subscribe(
        "/frankapy/gripper_grasp", 20, &CartesianPythonController::gripperGraspCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());
    
    this->_sub_joint = control_py_node.subscribe(
        "/frankapy/jointmotion", 20, &CartesianPythonController::jointMotionCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());

    this->_sub_joint_traj = control_py_node.subscribe(
        "/frankapy/joint_trajectory", 20, &CartesianPythonController::jointTrajectoryCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());
    

    this->_sub_joint_vel_traj = control_py_node.subscribe(
        "/frankapy/joint_vel_trajectory", 20, &CartesianPythonController::jointVelocityTrajectoryCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());

    this->_sub_move_relative = control_py_node.subscribe(
        "/frankapy/move_relative", 20, &CartesianPythonController::moveRelativeCallback, this,
        ros::TransportHints().reliable().tcpNoDelay());

    
    // set collision behavior
    this->_robot->setCollisionBehavior(
        {{100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0}},
        {{100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0}},
        {{100.0, 100.0, 100.0, 100.0, 100.0, 100.0}},
        {{100.0, 100.0, 100.0, 100.0, 100.0, 100.0}});
    this->_robot->setJointImpedance({{3000, 3000, 3000, 2500, 2500, 2000, 2000}});
    //this->_robot->setJointImpedance({{1000, 1000, 1000, 1000, 1000, 1000, 1000}});
    this->_robot->setCartesianImpedance({{3000, 3000, 3000, 300, 300, 300}});

}

void CartesianPythonController::jointMotionCallback(const franka_david::JointMotionPyPtr& msg)
{
    double joint0 = msg->joint0;
    if (this->_franka3)
    {
        joint0 = -msg->joint0;
    }

    double joint1 = msg->joint1;
    double joint2 = msg->joint2;
    if (this->_franka3)
    {
        joint2 = -msg->joint2;
    }
    double joint3 = msg->joint3;
    double joint4 = msg->joint4;
    if (this->_franka3)
    {
        joint4 = -msg->joint4;
    }
    double joint5 = msg->joint5;
    double joint6 = msg->joint6;
    if (this->_franka3)
    {
        joint6 = -msg->joint6 + M_PI/2; //old initial ee orientation
        //joint6 = msg->joint6; //new initial ee orientation
    }
    bool enable = msg->enable;
        
    if (enable)
    {
        std::array<double, 7> q_goal = {{joint0, joint1, joint2, joint3, joint4, joint5, joint6}};
        MotionGenerator motion_generator(0.2, q_goal);
        this->_robot->control(motion_generator);
    }

    
    
}

void CartesianPythonController::jointTrajectoryCallback(const franka_david::JointTrajPyPtr& msg)
{
    std::vector<double> joint0 =  msg->joint0; //list of joint0 values
    std::vector<double> joint1 =  msg->joint1; //list of joint1 values
    std::vector<double> joint2 =  msg->joint2; //list of joint2 values
    std::vector<double> joint3 =  msg->joint3; //list of joint3 values
    std::vector<double> joint4 =  msg->joint4; //list of joint4 values
    std::vector<double> joint5 =  msg->joint5; //list of joint5 values
    std::vector<double> joint6 =  msg->joint6; //list of joint6 values
    bool enable = msg->enable;
    double duration = msg->duration;

    //FLIP JOINTS 0, 2, 4 and 6 SIGNS and modify joint7 FOR FRANKA3
    if (this->_franka3)
    {
        std::transform(joint0.cbegin(),joint0.cend(),joint0.begin(),std::negate<double>());
        std::transform(joint2.cbegin(),joint2.cend(),joint2.begin(),std::negate<double>());
        std::transform(joint4.cbegin(),joint4.cend(),joint4.begin(),std::negate<double>());
        
        std::transform(joint6.cbegin(),joint6.cend(),joint6.begin(),std::negate<double>());
        std::transform(joint6.cbegin(),joint6.cend(),joint6.begin(),std::bind2nd(std::plus<double>(), M_PI/2));
    }
    
    std::cout << "Trajectory duration " << duration << std::endl;
    int N = int(duration / this->_dt);
    size_t loop_iter = 0;
    double time = 0.0;

    //std::cout << "iter: " << loop_iter << " sent joint values:"<< std::endl;
    //std::cout << "j1:" << joint0[loop_iter] << " j2:" << joint1[loop_iter] << " j3:" << joint2[loop_iter] << " j4:" << joint3[loop_iter] << " j5:" << joint4[loop_iter] << " j6:" << joint5[loop_iter] << " j7:" << joint6[loop_iter] <<  std::endl;

    if (enable)
    {
        //BASED ON https://petercorke.com/robotics/franka-emika-control-interface-libfranka/
        // and https://frankaemika.github.io/libfranka/joint_impedance_control_8cpp-example.html


        /*
        std::function<franka::JointPositions(const franka::RobotState &, franka::Duration)>
            joint_motion_generator_callback = [&](const franka::RobotState &robot_state,
                                             franka::Duration period) -> franka::JointPositions
        {
            


            std::cout << "iter: " << loop_iter << " sent joint values:"<< std::endl;
            std::cout << "j1:" << joint0[loop_iter] << " j2:" << joint1[loop_iter] << " j3:" << joint2[loop_iter] << " j4:" << joint3[loop_iter] << " j5:" << joint4[loop_iter] << " j6:" << joint5[loop_iter] << " j7:" << joint6[loop_iter] <<  std::endl;

            franka::JointPositions output = {{joint0[loop_iter], joint1[loop_iter], joint2[loop_iter],
            joint3[loop_iter], joint4[loop_iter], joint5[loop_iter], joint6[loop_iter]}};

            //ADDED TO PRINT FOLLOWED TRAJ TO FILE
            //get actual joint values as in https://frankaemika.github.io/libfranka/cartesian_impedance_control_8cpp-example.html#a12
            Eigen::Map<const Eigen::Matrix<double, 7, 1>> actual_q(robot_state.q.data());

            std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
            std::ofstream out_file(HOME + "/catkin_ws/src/Data/executed_joint_trajectory.csv", ios::app);
            out_file << actual_q[0] << "," << actual_q[1] << "," << actual_q[2] << "," << actual_q[3] << "," << actual_q[4]<< "," << actual_q[5] << "," << actual_q[6] << "," << std::endl;

            //loop_iter++;
            time += period.toSec();
            loop_iter = int(time / this->_dt);

            if(loop_iter > N-1){
                std::cout << std::endl << "Finished motion" << std::endl;
                return franka::MotionFinished(output);
            }

            return output;

        };

        this->_robot->control(joint_motion_generator_callback);
        //this->_robot->control(joint_motion_generator_callback, franka::ControllerMode::kJointImpedance);

        */

        

        std::array<double, 7> initial_position;
        this->_robot->control([&](const franka::RobotState& robot_state,
                                                franka::Duration period) -> franka::JointPositions {


        if (time == 0.0) {
            initial_position = robot_state.q_d;
        }
        double delta_angle = M_PI / 8.0 * (1 - std::cos(M_PI / 2.5 * time));
        franka::JointPositions output2 = {{initial_position[0], initial_position[1],
                                            initial_position[2], initial_position[3] + delta_angle,
                                            initial_position[4] + delta_angle, initial_position[5],
                                            initial_position[6] + delta_angle}};


        
        //std::cout << "j1:" << joint0[loop_iter] << " j2:" << joint1[loop_iter] << " j3:" << joint2[loop_iter] << " j4:" << joint3[loop_iter] << " j5:" << joint4[loop_iter] << " j6:" << joint5[loop_iter] << " j7:" << joint6[loop_iter] <<  std::endl;
        franka::JointPositions output = {{joint0[loop_iter], joint1[loop_iter], joint2[loop_iter],
            joint3[loop_iter], joint4[loop_iter], joint5[loop_iter], joint6[loop_iter]}};

        std::cout << "iter: " << loop_iter << std::endl;

        double diff0, diff1, diff2, diff3, diff4, diff5, diff6, diff7;
        diff0 = std::abs(initial_position[0] - joint0[loop_iter]);
        diff1 = std::abs(initial_position[1] - joint1[loop_iter]);
        diff2 = std::abs(initial_position[2] - joint2[loop_iter]);
        diff3 = std::abs(initial_position[3] - joint3[loop_iter]+ delta_angle);
        diff4 = std::abs(initial_position[4] - joint4[loop_iter]+ delta_angle);
        diff5 = std::abs(initial_position[5] - joint5[loop_iter]);
        diff6 = std::abs(initial_position[6] - joint6[loop_iter]+ delta_angle);

        double error_tolerance = 0.001;

        if ((diff0 > error_tolerance) || (diff1 > error_tolerance) || (diff2 > error_tolerance) || (diff3 > error_tolerance) || (diff4 > error_tolerance) || (diff5 > error_tolerance) || (diff6 > error_tolerance))
        {
        std::cout << "iter: " << loop_iter << " difference:"<< std::endl;
        std::cout << "j1:" <<initial_position[0] - joint0[loop_iter]   << " j2:" <<initial_position[1] - joint1[loop_iter] << " j3:" <<initial_position[2] - joint2[loop_iter];
        std::cout << " j4:" <<initial_position[3] - joint3[loop_iter] + delta_angle << " j5:" <<initial_position[4] - joint4[loop_iter] + delta_angle;
        std::cout << " j6:" <<initial_position[5] - joint5[loop_iter] << " j7:" <<initial_position[6] - joint6[loop_iter] + delta_angle << std::endl;
        }
        if (loop_iter > 1000 && loop_iter < 1100)
        {
            std::cout << "iter: " << loop_iter << " sent joint values:"<< std::endl;
            std::cout << "0 " << diff0 << " 1 " << diff1 << " 2 " << diff2 << " 3 " << diff3 << " 4 " << diff4 << " 5 " << diff5 << " 6 " << diff6 << std::endl;
        }


        //ADDED TO PRINT FOLLOWED TRAJ TO FILE
        //get actual joint values as in https://frankaemika.github.io/libfranka/cartesian_impedance_control_8cpp-example.html#a12
        Eigen::Map<const Eigen::Matrix<double, 7, 1>> actual_q(robot_state.q.data());

        std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
        std::ofstream out_file(HOME + "/catkin_ws/src/Data/executed_joint_trajectory.csv", ios::app);
        out_file << actual_q[0] << "," << actual_q[1] << "," << actual_q[2] << "," << actual_q[3] << "," << actual_q[4]<< "," << actual_q[5] << "," << actual_q[6] << "," << std::endl;
        //out_file << joint0[loop_iter] << "," << joint1[loop_iter] << "," << joint2[loop_iter] << "," << joint3[loop_iter] << "," << joint4[loop_iter] << "," << joint5[loop_iter] << "," << joint6[loop_iter] << std::endl;
        

        //joint velocities
        Eigen::Map<const Eigen::Matrix<double, 7, 1>> actual_dq(robot_state.dq.data());

        std::ofstream vel_out_file(HOME + "/catkin_ws/src/Data/executed_joint_velocities.csv", ios::app);
        vel_out_file << actual_dq[0] << "," << actual_dq[1] << "," << actual_dq[2] << "," << actual_dq[3] << "," << actual_dq[4]<< "," << actual_dq[5] << "," << actual_dq[6] << "," << std::endl;

        //Cartesian pose during joint control
        Eigen::Affine3d transform(Eigen::Matrix4d::Map(robot_state.O_T_EE.data()));
        Eigen::Vector3d position(transform.translation());
        Eigen::Quaterniond orientation(transform.linear()); // transform.linear() extracts the rotation matrix

        std::ofstream pose_out_file(HOME + "/catkin_ws/src/Data/joint_control_pose.csv", ios::app);
        pose_out_file << position[0] << "," << position[1] << "," << position[2] << "," << orientation.x() << "," << orientation.y() << "," << orientation.z() << "," << orientation.w() << "," << std::endl;



        time += period.toSec();
        std::cout << "Time "<< double(time) << std::endl;
        //loop_iter = int(time / this->_dt) + 1;
        loop_iter++;

        if (loop_iter > N-1) {
            std::cout << std::endl << "Finished motion, shutting down example" << std::endl;
            return franka::MotionFinished(output);
        }
        return output;
        });

    }
}

void CartesianPythonController::jointVelocityTrajectoryCallback(const franka_david::JointVelTrajPyPtr& msg)
{
    std::vector<double> joint0 =  msg->joint0; //list of joint0 values
    std::vector<double> joint1 =  msg->joint1; //list of joint1 values
    std::vector<double> joint2 =  msg->joint2; //list of joint2 values
    std::vector<double> joint3 =  msg->joint3; //list of joint3 values
    std::vector<double> joint4 =  msg->joint4; //list of joint4 values
    std::vector<double> joint5 =  msg->joint5; //list of joint5 values
    std::vector<double> joint6 =  msg->joint6; //list of joint6 values
    bool enable = msg->enable;
    double duration = msg->duration;

    //FLIP JOINTS 0, 2, 4 and 6 SIGNS and modify joint7 FOR FRANKA3
    if (this->_franka3)
    {
        std::transform(joint0.cbegin(),joint0.cend(),joint0.begin(),std::negate<double>());
        std::transform(joint2.cbegin(),joint2.cend(),joint2.begin(),std::negate<double>());
        std::transform(joint4.cbegin(),joint4.cend(),joint4.begin(),std::negate<double>());
        
        std::transform(joint6.cbegin(),joint6.cend(),joint6.begin(),std::negate<double>());
        //std::transform(joint6.cbegin(),joint6.cend(),joint6.begin(),std::bind2nd(std::plus<double>(), M_PI/2));
    }
    
    std::cout << "Trajectory duration " << duration << std::endl;
    int N = int(duration / this->_dt);
    size_t loop_iter = 0;
    double time = 0.0;

    if (enable)
    {
        //BASED ON https://petercorke.com/robotics/franka-emika-control-interface-libfranka/
        // and https://frankaemika.github.io/libfranka/joint_impedance_control_8cpp-example.html
        //FOR VELOCITY: https://frankaemika.github.io/libfranka/generate_joint_velocity_motion_8cpp-example.html

        this->_robot->control([&](const franka::RobotState& robot_state,
                                                franka::Duration period) -> franka::JointVelocities {


        //std::cout << "j1:" << joint0[loop_iter] << " j2:" << joint1[loop_iter] << " j3:" << joint2[loop_iter] << " j4:" << joint3[loop_iter] << " j5:" << joint4[loop_iter] << " j6:" << joint5[loop_iter] << " j7:" << joint6[loop_iter] <<  std::endl;
        franka::JointVelocities output = {{joint0[loop_iter], joint1[loop_iter], joint2[loop_iter],
            joint3[loop_iter], joint4[loop_iter], joint5[loop_iter], joint6[loop_iter]}};

        std::cout << "iter: " << loop_iter << std::endl;

        //ADDED TO PRINT FOLLOWED TRAJ TO FILE
        //get actual joint values as in https://frankaemika.github.io/libfranka/cartesian_impedance_control_8cpp-example.html#a12
        Eigen::Map<const Eigen::Matrix<double, 7, 1>> actual_q(robot_state.q.data());

        std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
        std::ofstream out_file(HOME + "/catkin_ws/src/Data/executed_joint_trajectory.csv", ios::app);
        out_file << actual_q[0] << "," << actual_q[1] << "," << actual_q[2] << "," << actual_q[3] << "," << actual_q[4]<< "," << actual_q[5] << "," << actual_q[6] << "," << std::endl;
        //out_file << joint0[loop_iter] << "," << joint1[loop_iter] << "," << joint2[loop_iter] << "," << joint3[loop_iter] << "," << joint4[loop_iter] << "," << joint5[loop_iter] << "," << joint6[loop_iter] << std::endl;
        

        //joint velocities
        Eigen::Map<const Eigen::Matrix<double, 7, 1>> actual_dq(robot_state.dq.data());

        std::ofstream vel_out_file(HOME + "/catkin_ws/src/Data/executed_joint_velocities.csv", ios::app);
        vel_out_file << actual_dq[0] << "," << actual_dq[1] << "," << actual_dq[2] << "," << actual_dq[3] << "," << actual_dq[4]<< "," << actual_dq[5] << "," << actual_dq[6] << "," << std::endl;

        //Cartesian pose during joint control
        Eigen::Affine3d transform(Eigen::Matrix4d::Map(robot_state.O_T_EE.data()));
        Eigen::Vector3d position(transform.translation());
        Eigen::Quaterniond orientation(transform.linear()); // transform.linear() extracts the rotation matrix

        std::ofstream pose_out_file(HOME + "/catkin_ws/src/Data/joint_control_pose.csv", ios::app);
        pose_out_file << position[0] << "," << position[1] << "," << position[2] << "," << orientation.x() << "," << orientation.y() << "," << orientation.z() << "," << orientation.w() << "," << std::endl;


        time += period.toSec();
        std::cout << "Time "<< double(time) << std::endl;
        //loop_iter = int(time / this->_dt) + 1;
        loop_iter++;

        if (loop_iter > N-1) {
            std::cout << std::endl << "Finished motion, shutting down example" << std::endl;
            return franka::MotionFinished(output);
        }
        return output;
        });

    }
}

void CartesianPythonController::gripperMoveCallback(const franka_david::GripperPyPtr& msg)
{
    double width = msg->width;
    double speed = msg->speed;
    bool enable = msg->enable;
    
    if (enable)
    {
    	std::cout << "Moving the gripper at speed " << speed << " and to width " << width << std::endl;
        this->_gripper->move(width, speed);
    }
    
}

void CartesianPythonController::gripperGraspCallback(const franka_david::GripperGraspPyPtr& msg)
{
    double distance = msg->distance;
    double speed = msg->speed;
    double force = msg->force;
    bool release = msg->release;
    bool enable = msg->enable;
    
    if (enable)
    {
    	
        if (release)
        {
            std::cout << "Releasing gripper" <<std::endl;
            this->_gripper->stop();
        }
        else
        {
            std::cout << "Grasping with force " << force << "N and distance" << distance << std::endl;
            this->_gripper->grasp(distance, speed, force);
        }
    }
    
}



void CartesianPythonController::movetoCallback(const franka_david::MovetoPyPtr& msg)
{   
    
    geometry_msgs::Point pose = msg->pose.pose.position;
    geometry_msgs::Quaternion ori_d = msg->pose.pose.orientation;
    double duration = msg->duration;
    double pose_x = pose.x;
    double pose_y = pose.y;
    double pose_z = pose.z;
    double ori_x = ori_d.x;
    double ori_y = ori_d.y;
    double ori_z = ori_d.z;
    double ori_w = ori_d.w;
    bool enable = msg->enable;
    
    std::cout << "Trajectory duration " << duration << std::endl;
    // time duration in seconds and steps
//     double duration = 4;
    int N = int(duration / this->_dt);
    int Nend = N * 10 /100;  // Percentage of the trajectory to determine point from which to have more trajectory points
    int Npoints1 = N * 60/100;  // Number of points designated for the first part of the traj
//     std::cout << "Trajectory length=" << N << ", end points=" << Nend << std::endl;
    
    Eigen::Affine3d transform;
    Eigen::Vector3d pos;
    Eigen::Quaterniond ori;
    math::Transform3D *Ti_temp, *Ti_temp1, *Ti_temp2, *Ti, T0, Tmiddle, Tend;
    // Read the current state of the robot
    
    if (enable)
    {
    	this->_robot->read([&](const franka::RobotState &robot_state){
		transform = Eigen::Matrix4d::Map(robot_state.O_T_EE.data());
		pos = transform.translation();
		ori = transform.linear();
		Eigen::Map<const Eigen::Matrix<double, 7, 1>> joint_angles(robot_state.q.data());
		Eigen::Matrix3d R = transform.rotation();
// 		T0(1, 1) = R(0, 0);    T0(1, 2) = R(0, 1);    T0(1, 3) = R(0, 2);    T0(1, 4) = pos[0];
// 		T0(2, 1) = R(1, 0);    T0(2, 2) = R(1, 1);    T0(2, 3) = R(1, 2);    T0(2, 4) = pos[1];
// 		T0(3, 1) = R(2, 0);    T0(3, 2) = R(2, 1);    T0(3, 3) = R(2, 2);    T0(3, 4) = pos[2];
        
//         T0(1, 1) = 0.9981;    T0(1, 2) = -0.0315;    T0(1, 3) = 0.0531;    T0(1, 4) = pos[0];
// 		T0(2, 1) = -0.0315;    T0(2, 2) = -0.9995;    T0(2, 3) = -0.0002;    T0(2, 4) = pos[1];
// 		T0(3, 1) = 0.0531;    T0(3, 2) = -0.0015;    T0(3, 3) = -0.9986;    T0(3, 4) = pos[2];
        
        
        T0(1, 1) = 1.0;    T0(1, 2) = -0.0;    T0(1, 3) = 0.0;    T0(1, 4) = pos[0];
		T0(2, 1) = -0.0;    T0(2, 2) = -1.0;    T0(2, 3) = -0.0;    T0(2, 4) = pos[1];
		T0(3, 1) = 0.0;    T0(3, 2) = -0.0;    T0(3, 3) = -1.0;    T0(3, 4) = pos[2];

		Tend = T0;
		Tend(1, 4) = pose_x;
		Tend(2, 4) = pose_y;
		Tend(3, 4) = pose_z;

		T0.print("my current");
		Tend.print("my goal");

        //ADDED TO PRINT FOLLOWED TRAJ TO FILE
        // std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
        // std::ofstream out_file(HOME + "/catkin_ws/src/franka_david/scripts/real_traj/output.csv", ios::app);    
        // out_file << pos[0] << "," << pos[1] << "," << pos[2] << "," << ori.x() << "," << ori.y() << "," << ori.z() << "," << ori.w() << "," << std::endl;


        // if (!out_file)
        // {
        //     cout << "Error in creating file!!!" << endl; //TURNS OUT THAT FILE CREATION FAILS > why?
        // }

        // out_file.close();

		return 0;
	    });
        
        
	    Ti_temp = math::ctraj(T0, Tend, N);
	    // Tmiddle = T0;
	    // Tmiddle(1, 4) = Ti_temp[N-Nend](1, 4);
	    // Tmiddle(2, 4) = Ti_temp[N-Nend](2, 4);
	    // Tmiddle(3, 4) = Ti_temp[N-Nend](3, 4);
	    // //std::cout << "Compute T temp 1" << std::endl;
	    
	    // Ti_temp1 = math::ctraj(T0, Tmiddle, Npoints1);
	    // //std::cout << "Compute T temp 2" << std::endl;
	    // Ti_temp2 = math::ctraj(Tmiddle, Tend, N-Npoints1);
	    // //std::cout << "For loop" << std::endl;
	    
	    // Ti = new math::Transform3D[(int)N];
	    
	    // for (int k=0; k<(int)N; k++)
	    // {   
	    //     if (k< Npoints1)
	    //     {
	    //         Ti[k] = Ti_temp1[k];
	    //         //std::cout << "iteration" << k << std::endl;
	    //     }
	    //     else
	    //     {
	    //         //std::cout << "Number of points" << k-(Npoints1) << std::endl;
	    //         Ti[k] = Ti_temp2[k-(Npoints1)];
	    //     }
	        
	    // }

	    
	    this->runControl(Ti_temp, N);
    }
    
    
}


void CartesianPythonController::motionCallback(const franka_david::MotionPyPtr& msg)
{
        
    std::vector<geometry_msgs::PoseStamped> trajectory = msg->trajectory;
    bool enable = msg->enable;
    
    if (enable)
    {   
        std::vector<geometry_msgs::PoseStamped>::iterator it;
    
        std::vector<double> trajectory_x;
        std::vector<double> trajectory_y;
        std::vector<double> trajectory_z;
        

        for(it = trajectory.begin(); it != trajectory.end(); it++)    {
            
            trajectory_x.push_back(it->pose.position.x);
            trajectory_y.push_back(it->pose.position.y);
            trajectory_z.push_back(it->pose.position.z);

        }
        int N = trajectory_x.size();
        std::cout << "trajectory size=" << N << std::endl;
        double duration = N * this->_dt;
        std::cout << "Trajectory duration " << duration << std::endl;
        double end_pose_x = trajectory_x.at(N-1);
        double end_pose_y = trajectory_y.at(N-1);
        double end_pose_z = trajectory_z.at(N-1);
        
        
        Eigen::Affine3d transform;
        Eigen::Vector3d pos;
        Eigen::Quaterniond ori;
        math::Transform3D *Ti, T0, Tend;

        
    	this->_robot->read([&](const franka::RobotState &robot_state){
		transform = Eigen::Matrix4d::Map(robot_state.O_T_EE.data());
		pos = transform.translation();
		ori = transform.linear();
		Eigen::Map<const Eigen::Matrix<double, 7, 1>> joint_angles(robot_state.q.data());
		Eigen::Matrix3d R = transform.rotation();

//         T0(1, 1) = 0.9981;    T0(1, 2) = -0.0315;    T0(1, 3) = 0.0531;    T0(1, 4) = pos[0];
// 		T0(2, 1) = -0.0315;    T0(2, 2) = -0.9995;    T0(2, 3) = -0.0002;    T0(2, 4) = pos[1];
// 		T0(3, 1) = 0.0531;    T0(3, 2) = -0.0015;    T0(3, 3) = -0.9986;    T0(3, 4) = pos[2];
        T0(1, 1) = 1.0;    T0(1, 2) = -0.0;    T0(1, 3) = 0.0;    T0(1, 4) = pos[0];
		T0(2, 1) = -0.0;    T0(2, 2) = -1.0;    T0(2, 3) = -0.0;    T0(2, 4) = pos[1];
		T0(3, 1) = 0.0;    T0(3, 2) = -0.0;    T0(3, 3) = -1.0;    T0(3, 4) = pos[2];

		Tend = T0;
		Tend(1, 4) = end_pose_x;
		Tend(2, 4) = end_pose_y;
		Tend(3, 4) = end_pose_z;

		T0.print("my current");
		Tend.print("my goal");


		return 0;
	    });
        
	    Ti = math::ctraj(T0, Tend, N);
	    
	    for (int k=0; k<(int)N; k++)
	    {   
            Ti[k](1, 4) = trajectory_x.at(k);
            Ti[k](2, 4) = trajectory_y.at(k);
            Ti[k](3, 4) = trajectory_z.at(k);
// 	        std::cout << "X=" << trajectory_x.at(k) << ", Y=" << trajectory_y.at(k) << ", Z=" << trajectory_z.at(k) << std::endl;
	    }
	    
	    
	    this->runControl(Ti, N);
    }
    
    
    
// //     geometry_msgs::Quaternion ori_d = msg->pose.pose.orientation;
//     double duration = msg->duration;
// 
//     bool enable = msg->enable;
//     
// 
//     
}


void CartesianPythonController::motionOriCallback(const franka_david::MotionPyPtr& msg)
{
        
    std::vector<geometry_msgs::PoseStamped> trajectory = msg->trajectory;
    bool enable = msg->enable;
    
    if (enable)
    {   
        std::vector<geometry_msgs::PoseStamped>::iterator it;
    
        std::vector<double> trajectory_x;
        std::vector<double> trajectory_y;
        std::vector<double> trajectory_z;
        std::vector<Eigen::Quaterniond> trajectory_quat;


        for(it = trajectory.begin(); it != trajectory.end(); it++)    {
            
            trajectory_x.push_back(it->pose.position.x);
            trajectory_y.push_back(it->pose.position.y);
            trajectory_z.push_back(it->pose.position.z);
            Eigen::Quaterniond a;
            a.x() = it->pose.orientation.x;
            a.y() = it->pose.orientation.y;
            a.z() = it->pose.orientation.z;
            a.w() = it->pose.orientation.w;
            trajectory_quat.push_back(a);

        }
        int N = trajectory_x.size();
        std::cout << "trajectory size=" << N << std::endl;
        double duration = N * this->_dt;
        std::cout << "Trajectory duration " << duration << std::endl;
        double end_pose_x = trajectory_x.at(N-1);
        double end_pose_y = trajectory_y.at(N-1);
        if (this->_franka3)
        {
            end_pose_y = -trajectory_y.at(N-1);
        }
        double end_pose_z = trajectory_z.at(N-1);
        
        
        Eigen::Affine3d transform;
        Eigen::Vector3d pos;
        Eigen::Quaterniond ori;
        math::Transform3D *Ti, T0, Tend;

        
    	this->_robot->read([&](const franka::RobotState &robot_state){
		transform = Eigen::Matrix4d::Map(robot_state.O_T_EE.data());
		pos = transform.translation();
		ori = transform.linear();
		// Eigen::Map<const Eigen::Matrix<double, 7, 1>> joint_angles(robot_state.q.data());
		Eigen::Matrix3d R = transform.rotation();

//         T0(1, 1) = 0.9981;    T0(1, 2) = -0.0315;    T0(1, 3) = 0.0531;    T0(1, 4) = pos[0];
// 		T0(2, 1) = -0.0315;    T0(2, 2) = -0.9995;    T0(2, 3) = -0.0002;    T0(2, 4) = pos[1];
// 		T0(3, 1) = 0.0531;    T0(3, 2) = -0.0015;    T0(3, 3) = -0.9986;    T0(3, 4) = pos[2];
        if (this->_franka3)
        {
            T0(1, 1) = 1.0;    T0(1, 2) = -0.0;    T0(1, 3) = 0.0;    T0(1, 4) = pos[0];
            T0(2, 1) = -0.0;    T0(2, 2) = -1.0;    T0(2, 3) = -0.0;    T0(2, 4) = pos[1];
            T0(3, 1) = 0.0;    T0(3, 2) = -0.0;    T0(3, 3) = -1.0;    T0(3, 4) = pos[2];
        }
        else
        {
            T0(1, 1) = 1.0;    T0(1, 2) = -0.0;    T0(1, 3) = 0.0;    T0(1, 4) = pos[0];
            T0(2, 1) = -0.0;    T0(2, 2) = -1.0;    T0(2, 3) = -0.0;    T0(2, 4) = pos[1];
            T0(3, 1) = 0.0;    T0(3, 2) = -0.0;    T0(3, 3) = -1.0;    T0(3, 4) = pos[2];
        }

		Tend = T0;
		Tend(1, 4) = end_pose_x;
		Tend(2, 4) = end_pose_y;
		Tend(3, 4) = end_pose_z;

		T0.print("my current");
		Tend.print("my goal");


		return 0;
	    });
        
	    Ti = math::ctraj(T0, Tend, N);
	    
	    for (int k=0; k<(int)N; k++)
	    {   
            Ti[k](1, 4) = trajectory_x.at(k);
            if (this->_franka3)
            {
                Ti[k](2, 4) = -trajectory_y.at(k);
            }
            else
            {
                Ti[k](2, 4) = trajectory_y.at(k);
            }
            Ti[k](3, 4) = trajectory_z.at(k);
// 	        std::cout << "X=" << trajectory_x.at(k) << ", Y=" << trajectory_y.at(k) << ", Z=" << trajectory_z.at(k) << std::endl;
	        Eigen::Matrix3d rot = trajectory_quat.at(k).normalized().toRotationMatrix();
            Ti[k](1, 1) = rot(0, 0); Ti[k](1, 2) = rot(0, 1); Ti[k](1, 3) = rot(0, 2); 
            Ti[k](2, 1) = rot(1, 0); Ti[k](2, 2) = rot(1, 1); Ti[k](2, 3) = rot(1, 2); 
            Ti[k](3, 1) = rot(2, 0); Ti[k](3, 2) = rot(2, 1); Ti[k](3, 3) = rot(2, 2);

            //if (! this->_franka3) //OLD
            if (this->_franka3)
            {
                //From OLD code for other inital ee rotation
                // Ti[k](1, 2) = -rot(0, 1);
                // Ti[k](2, 1) = -rot(1, 0);
                // Ti[k](2, 3) = -rot(1, 2); 
                // Ti[k](3, 2) = -rot(2, 1);

                //For NEW ee initial rotation - NOTE 04.09: also seems to work for original ee orientation
                Ti[k](3, 1) = -rot(2, 0);
                Ti[k](3, 2) = -rot(2, 1);
                Ti[k](1, 3) = -rot(0, 2); 
                Ti[k](2, 3) = -rot(1, 2);
            }
        }
	    
	    this->runControl(Ti, N);
    }
    
    
    
// //     geometry_msgs::Quaternion ori_d = msg->pose.pose.orientation;
//     double duration = msg->duration;
// 
//     bool enable = msg->enable;
//     
// 
//     
}


void CartesianPythonController::moveRelativeCallback(const franka_david::MoveRelativePyPtr& msg)
{
    double dx = msg->x;
    double dy = msg->y;
    if (this->_franka3)
    {
        dy = -dy;
    }
    double dz = msg->z;
    double duration = msg->duration;
    bool enable = msg->enable;
        
    
    if (enable)
    {
    //Motion generator based on this example: https://frankaemika.github.io/libfranka/generate_cartesian_pose_motion_8cpp-example.html
    std::array<double, 16> initial_pose;
    double time = 0.0;

    this->_robot->control([&](const franka::RobotState& robot_state,
                                         franka::Duration period) -> franka::CartesianPose {
        if (time == 0.0) {
        initial_pose = robot_state.O_T_EE_c;
        }
        time += period.toSec();
        //double delta_x = (time/duration)*dx; //Old, now instead of linear velocity use the lines below to slow down velocity at end of motion
        double angle = M_PI / 4 * (1 - std::cos(M_PI / duration * time));
        double delta_x = dx * std::sin(angle);
        double delta_y = dy * std::sin(angle);
        double delta_z = dz * std::sin(angle);
        
        std::array<double, 16> new_pose = initial_pose;
        new_pose[12] += delta_x;
        new_pose[13] += delta_y;
        new_pose[14] += delta_z;


        if (time >= duration) {
            
            std::cout << std::endl << "Finished motion, shutting down example" << std::endl;
            return franka::MotionFinished(new_pose);
        }
        return new_pose;
    });
    }
}


void CartesianPythonController::runControl(math::Transform3D* trajectory, int N)
{
    // Compliance parameters
    // 2000 - smoother than 3000
    const double translational_stiffness{1400.0}; //original 2000. 150 works nice for 1000fps Franka3. 108 Franka2.
//     100 shakes too much, 50 not that good, 10- too soft
    // 70 might be fine - shakes a bit, 
    // 68 shakes close to the robot, doesn't when far but orientation is really good
    const double rotational_stiffness{53.0}; //original 68. 40 works nice for 1000fps Franka3 and Franka2.
    Eigen::MatrixXd stiffness(6, 6), damping(6, 6);
    stiffness.setZero();
    stiffness.topLeftCorner(3, 3) << translational_stiffness * Eigen::MatrixXd::Identity(3, 3);
    stiffness.bottomRightCorner(3, 3) << rotational_stiffness * Eigen::MatrixXd::Identity(3, 3);
    damping.setZero();
    damping.topLeftCorner(3, 3) << 2.0 * sqrt(translational_stiffness) *
                                    Eigen::MatrixXd::Identity(3, 3);
    damping.bottomRightCorner(3, 3) << 2.0 * sqrt(rotational_stiffness) *
                                        Eigen::MatrixXd::Identity(3, 3);
                                        
    size_t loop_iter = 0;
    
    Eigen::Matrix<double, 6, 1> error, error_old;
    error_old << .0, .0, .0, .0, .0, .0;
    Eigen::Vector3d position_d;
    Eigen::Quaterniond orientation_d;
    franka::Model model = this->_robot->loadModel();
    
    double dt = this->_dt;
    double time = 0.0;
    
    std::function<franka::Torques(const franka::RobotState &, franka::Duration)>
            impedance_control_callback = [&](const franka::RobotState &robot_state,
                                             franka::Duration /*duration*/) -> franka::Torques
        {
            // get state variables
            std::array<double, 7> coriolis_array = model.coriolis(robot_state);
            std::array<double, 42> jacobian_array =
                model.zeroJacobian(franka::Frame::kEndEffector, robot_state);

            // convert to Eigen
            Eigen::Map<const Eigen::Matrix<double, 7, 1>> coriolis(coriolis_array.data());
            Eigen::Map<const Eigen::Matrix<double, 6, 7>> jacobian(jacobian_array.data());
            Eigen::Map<const Eigen::Matrix<double, 7, 1>> q(robot_state.q.data());
            Eigen::Map<const Eigen::Matrix<double, 7, 1>> dq(robot_state.dq.data());
            Eigen::Map<const Eigen::Matrix<double, 7, 1>> tau_J_d(robot_state.tau_J_d.data());
            Eigen::Affine3d transform(Eigen::Matrix4d::Map(robot_state.O_T_EE.data()));
            Eigen::Vector3d position(transform.translation());
            Eigen::Quaterniond orientation(transform.linear());
            // Desired trajectory
            math::Quaternion quat_desired(math::T2Q(trajectory[loop_iter]));
            orientation_d.w() = quat_desired.s();
            orientation_d.x() = quat_desired.x();
            orientation_d.y() = quat_desired.y();
            orientation_d.z() = quat_desired.z();		
            position_d << trajectory[loop_iter](1,4), trajectory[loop_iter](2,4), trajectory[loop_iter](3,4);
            // compute error to desired equilibrium pose
            // position error
            Eigen::Matrix<double, 6, 1> error;
            error.head(3) << position - position_d;
            // Prevent quaternion flips
            if (orientation_d.coeffs().dot(orientation.coeffs()) < 0.0) {
                orientation.coeffs() << -orientation.coeffs();
            }	
            // "difference" quaternion
            Eigen::Quaterniond error_quaternion(orientation.inverse() * orientation_d);
            error.tail(3) << error_quaternion.x(), error_quaternion.y(), error_quaternion.z();
            // Transform to base frame
            error.tail(3) << -transform.linear() * error.tail(3);
            // compute control
            Eigen::VectorXd tau_task(7), tau_d(7), force_task(6);

            // Spring damper system with damping ratio=1
            force_task << (-stiffness * error - damping * ((error - error_old)/dt));
            force_task << (-stiffness * error - damping * (jacobian * dq));
            tau_task << jacobian.transpose() * force_task;
            error_old = error;
            tau_d << tau_task + coriolis;
            // Saturate torque rate to avoid discontinuities
            tau_d = saturateTorqueRate(tau_d, tau_J_d);

            franka::Torques tau_d_array = {{tau_d(0), tau_d(1), tau_d(2), tau_d(3), tau_d(4), tau_d(5), tau_d(6)}};
            //std::array<double, 7> tau_d_array{};
            //Eigen::VectorXd::Map(&tau_d_array[0], 7) = tau_d;

            //ADDED TO PRINT FOLLOWED TRAJ TO FILE
            std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
            std::ofstream out_file(HOME + "/catkin_ws/src/Data/executed_trajectory.csv", ios::app);
            // if (!out_file)
            // {
            //     cout << "Error in creating file!!!" << endl; //TURNS OUT THAT FILE CREATION FAILS > why?
            // }
            out_file << position[0] << "," << position[1] << "," << position[2] << "," << orientation.x() << "," << orientation.y() << "," << orientation.z() << "," << orientation.w() << "," << std::endl;


            Eigen::Map<const Eigen::Matrix<double, 7, 1>> actual_q(robot_state.q.data());

            std::ofstream out_file_joints(HOME + "/catkin_ws/src/Data/pose_control_joints.csv", ios::app);
            out_file_joints << actual_q[0] << "," << actual_q[1] << "," << actual_q[2] << "," << actual_q[3] << "," << actual_q[4]<< "," << actual_q[5] << "," << actual_q[6] << "," << std::endl;


            loop_iter++;
            time += dt;
            if(loop_iter > N - 1) {
                std::cout << std::endl << "Finished motion" << std::endl;
                return franka::MotionFinished(tau_d_array);
            }
            return tau_d_array;
        };
        
    // start real-time control loop
//     std::cout << "WARNING: Collision thresholds are set to high values. "
//                 << "Make sure you have the user stop at hand!" << std::endl
//                 << "After starting try to push the robot and see how it reacts." << std::endl
//                 << "Press Enter to continue..." << std::endl;
//     std::cin.ignore();
    this->_robot->control(impedance_control_callback);
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "py_control");
    CartesianPythonController controller;

    ros::spin();

    return 0;

}

