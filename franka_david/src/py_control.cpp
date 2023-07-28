// Copyright (c) 2017 Franka Emika GmbH
// Use of this source code is governed by the Apache-2.0 license, see LICENSE
// #include <ros/ros.h>
// #include <ros/package.h>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>

#include <Eigen/Dense>

#include <franka/duration.h>
#include <franka/exception.h>
#include <franka/model.h>
#include <franka/robot.h>

#include "matrix.hpp"
#include "Frame3D.h"
#include "quaternion.hpp"

#include <examples_common.h>

#include <franka/exception.h>
#include <franka/gripper.h>


/**
 * @example cartesian_impedance_control.cpp
 * An example showing a simple cartesian impedance controller without inertia shaping
 * that renders a spring damper system where the equilibrium is the initial configuration.
 * After starting the controller try to push the robot around and try different stiffness levels.
 *
 * @warning collision thresholds are set to high values. Make sure you have the user stop at hand!
 */


int main(int argc, char** argv) {
    // Check whether the required arguments were passed

    std::string robot_hostname = "172.168.0.2";
	bool run = true;
	std::string c; 
		
    Eigen::Affine3d transform;
    Eigen::Vector3d pos;
    Eigen::Quaterniond ori;
    math::Transform3D *Ti, T0, T1;

    std::cin >> c;
    try {
        franka::Robot robot(robot_hostname);
        setDefaultBehavior(robot);

        // set collision behavior
        robot.setCollisionBehavior({{100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0}},
                                    {{100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0}},
                                    {{100.0, 100.0, 100.0, 100.0, 100.0, 100.0}},
                                    {{100.0, 100.0, 100.0, 100.0, 100.0, 100.0}});

        std::cout << "WARNING: This example will move the robot! "
                    << "Please make sure to have the user stop button at hand!" << std::endl
                    << "Press Enter to continue..." << std::endl;
        std::cin.ignore();

        // time duration and steps
        double duration = 5, dt = 0.001;
        int N = int(duration / dt);
        double time = 0.0;
        // define callback for the torque control loop
        size_t loop_iter = 0;

        // Load demonstrated trajectory
        int pose_dim = 7;
        std::vector<std::vector<double>> pqPath, stif_data, damp_data;
        std::string packPath = ros::package::getPath("franka_fares");
        int fileFound = loadVectorMatrixFromFile(packPath+"/in/spl_pos_uq.txt", pose_dim, pqPath);
        std::cout << "pqPath size: (" << pqPath.size() << "," << pqPath[0].size() << ")" << std::endl;
        size_t number_of_points = pqPath.size();

        // Allocate space to store execution data
        std::vector< std::vector<double> > execution, flow_traj;
        for(size_t i=0; i<number_of_points; ++i) {
            execution.push_back(std::vector<double>());
            for(size_t j=0; j<26; ++j) {
                execution[i].push_back(0.0);
            }
        }

        // load the kinematics and dynamics model
        franka::Model model = robot.loadModel();

        // equilibrium point is the initial position
        franka::RobotState initial_state = robot.readOnce();
        Eigen::Affine3d initial_transform(Eigen::Matrix4d::Map(initial_state.O_T_EE.data()));
        Eigen::Vector3d position_d(initial_transform.translation());
        Eigen::Quaterniond orientation_d(initial_transform.linear());
        double init_quat[4] = {orientation_d.w(), orientation_d.x(), orientation_d.y(), orientation_d.z()};
        // set nullspace equilibrium configuration to initial q
        Eigen::Map<Eigen::Matrix<double, 7, 1>> q_d_nullspace_(initial_state.q.data());
        
        // Compliance parameters
        const double translational_stiffness{2000.0};
        const double rotational_stiffness{70.0};
        Eigen::MatrixXd stiffness(6, 6), damping(6, 6);
        stiffness.setZero();
        stiffness.topLeftCorner(3, 3) << translational_stiffness * Eigen::MatrixXd::Identity(3, 3);
        stiffness.bottomRightCorner(3, 3) << rotational_stiffness * Eigen::MatrixXd::Identity(3, 3);
        damping.setZero();
        damping.topLeftCorner(3, 3) << 2.0 * sqrt(translational_stiffness) *
                                        Eigen::MatrixXd::Identity(3, 3);
        damping.bottomRightCorner(3, 3) << 2.0 * sqrt(rotational_stiffness) *
                                            Eigen::MatrixXd::Identity(3, 3);
        
        Eigen::Matrix<double, 6, 1> error, error_old;
        error_old << .0, .0, .0, .0, .0, .0;
        // define callback for the torque control loop
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
            orientation_d.w() = pqPath[loop_iter][3];
            orientation_d.x() = pqPath[loop_iter][4];
            orientation_d.y() = pqPath[loop_iter][5];
            orientation_d.z() = pqPath[loop_iter][6];		
            position_d << pqPath[loop_iter][0], pqPath[loop_iter][1], pqPath[loop_iter][2];

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

            // Save executed and generated trajectories	
            execution[loop_iter][0] = position[0];
            execution[loop_iter][1] = position[1];
            execution[loop_iter][2] = position[2];
            execution[loop_iter][3] = orientation.w();
            execution[loop_iter][4] = orientation.x();
            execution[loop_iter][5] = orientation.y();
            execution[loop_iter][6] = orientation.z();
            // External Wrenches
            execution[loop_iter][7] = robot_state.O_F_ext_hat_K[0];
            execution[loop_iter][8] = robot_state.O_F_ext_hat_K[1];
            execution[loop_iter][9] = robot_state.O_F_ext_hat_K[2];
            execution[loop_iter][10] = robot_state.O_F_ext_hat_K[3];
            execution[loop_iter][11] = robot_state.O_F_ext_hat_K[4];
            execution[loop_iter][12] = robot_state.O_F_ext_hat_K[5];
            // Control Wrenches
            execution[loop_iter][13] = force_task(0);
            execution[loop_iter][14] = force_task(1);
            execution[loop_iter][15] = force_task(2);
            execution[loop_iter][16] = force_task(3);
            execution[loop_iter][17] = force_task(4);
            execution[loop_iter][18] = force_task(5);
            // Control torques
            execution[loop_iter][19] = tau_task(0);
            execution[loop_iter][20] = tau_task(1);
            execution[loop_iter][21] = tau_task(2);
            execution[loop_iter][22] = tau_task(3);
            execution[loop_iter][23] = tau_task(4);
            execution[loop_iter][24] = tau_task(5);
            execution[loop_iter][25] = tau_task(6);
            // Joint angles
            execution[loop_iter][26] = q(0);
            execution[loop_iter][27] = q(1);
            execution[loop_iter][28] = q(2);
            execution[loop_iter][29] = q(3);
            execution[loop_iter][30] = q(4);
            execution[loop_iter][31] = q(5);
            execution[loop_iter][32] = q(6);
                
            loop_iter++;
            time += dt;
            if(loop_iter > number_of_points - 1) {
                std::cout << std::endl << "Finished motion, shutting down example" << std::endl;
                // release object.
                franka::Gripper gripper(argv[1]);
                gripper.stop();
                saveVectorMatrixToFile(packPath+"/out/riemannflow.txt", execution);
                return franka::MotionFinished(tau_d_array);
            }
            return tau_d_array;
        };
        // start real-time control loop
        std::cout << "WARNING: Collision thresholds are set to high values. "
                << "Make sure you have the user stop at hand!" << std::endl
                << "After starting try to push the robot and see how it reacts." << std::endl
                << "Press Enter to continue..." << std::endl;
        std::cin.ignore();
        robot.control(impedance_control_callback);
        // Save execution to file
        std::cout << std::endl << "Execution finished. Do you want to save it to file? [y/n]" << std::endl;
        string user_in;
        std::cin >> user_in;
        if(user_in=="y" || user_in=="Y") {
            std::cout << "Provide the execution number" << std::endl;
            std::string demoNum;
            std::cin >> demoNum;
            saveVectorMatrixToFile(packPath+"/out/riemannflow_"+demoNum+".txt", execution);
        }
    } catch (const franka::Exception& e) {
        std::cout << e.what() << std::endl;
        run = false;  
    }

  return 0;
}
