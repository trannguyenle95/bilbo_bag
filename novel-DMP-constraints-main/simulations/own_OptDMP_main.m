%NOTE: modifications also to other files, e.g. in offlineGMPweigths number
%of DMP kernels are no longer reset to 30

clc;
close all;
clear;

addpath('utils/');
import_gmp_lib();
import_io_lib();

addpath('/home/erichannus/catkin_ws/src/SupportScripts/')
addpath('/home/erichannus/catkin_ws/src/Data/demos/')

%% Load training data
bag = 'B';
filename = '10l_bag_flip.csv';

if strcmp('A', bag)
    bagwidth = 0.44
elseif strcmp('B', bag)
    bagwidth = 0.37
elseif strcmp('C', bag)
    bagwidth = 0.55
elseif strcmp('D', bag)
    bagwidth = 0.49
elseif strcmp('E', bag)
    bagwidth = 0.40
end


version = 4 % 3 = vel, 4 = pos, select which optimization version to export

D = preprocess(filename, false, 0.00, 0.00, 0.00, 1, 'ori1', bagwidth);
Dsmooth = smoothdata(D, 1, "gaussian", 35); %smooth demo before calculating IK
Dsmooth(:,4:7) = Dsmooth(:,4:7) ./ sqrt(sum(Dsmooth(:,4:7).^2,2)); %Make sure quaternion still has unit norm
[q, jacobians] = InverseKinematics(Dsmooth);
demo_traj = generateDemo(q', 1/120);

Timed = demo_traj.t;
Pd_data = demo_traj.pos;
dPd_data = demo_traj.vel;
ddPd_data = demo_traj.acc;

%% initialize and train GMP
train_method = 'LS';
N_kernels = 40; %Originally 30, difference in MSE but not noticable in plots
kernels_std_scaling = 1.5;
n_dof = size(Pd_data,1);
gmp = GMP(n_dof, N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
offline_train_mse
toc

gmp.setScaleMethod( TrajScale_Prop(n_dof) );

taud = Timed(end);
yd0 = gmp.getYd(0); %Pd_data(1);
ygd = gmp.getYd(1); %Pd_data(end);

kt = 1; %no temporal scaling of demo traj
%ks = diag([1 1 1]); % spatial scaling
tau = taud/kt;
y0 = yd0; %intial state same as as in demo
yg = ygd; %rollout goal is same as demo goal


%% ======== Limits ==========

%limits are chosen to be worst case between Franka Panda and Franka
%Research 3
pos_lim = [[-2.7437; -1.7628; -2.8973; -3.0421; -2.8065; 0.5445; -2.8973] [2.7437; 1.7628; 2.8973; -0.1518; 2.8065; 3.7525; 2.8973]];
vel_lim = [-[2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100] [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100]];
accel_lim = [-[10; 7.5; 10; 10; 10; 10; 10] [10; 7.5; 10; 10; 10; 10; 10]];

%use the actual robot limits for plotting
actual_pos_lim = pos_lim;
actual_vel_lim = vel_lim;
actual_accel_lim = accel_lim;

%use stricter limits than the actual ones for trajectory generation to have
%some margin
pos_lim = 1.0 * pos_lim; %no margin for position limit, so it is similar to TC-DMP which doesn't constrain pos but instead just fits DMP to position demo
vel_lim = 0.95 * vel_lim;
accel_lim = 0.95 * accel_lim;

%% ======== Generate trajectories ==========

qp_solver_type = 1; % matlab-quadprog:0 , osqp:1, Goldfarb-Idnani: 2
data = {};

%% --------- Demo -----------
data{length(data)+1} = ...
    struct('Time',Timed, 'Pos',Pd_data, 'Vel',dPd_data, 'Accel',ddPd_data, 'linestyle','--', ...
    'color',[0.7 0.7 0.7], 'legend','demo', 'plot3D',true, 'plot2D',true);

%% --------- Unconstrained DMP -----------
gmp2 = gmp.deepCopy();
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp2, tau, y0, yg);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
    'color','blue', 'legend','DMP', 'plot3D',true, 'plot2D',true);

%output unconstrained DMP for parameter comparison
writematrix(data{2}.Pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat(bag,'_UC_Opt_joint_',filename)))
writematrix(data{2}.Vel',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat(bag,'_UC_Opt_joint_vel_',filename)))

%% --------- Optimized DMP -> VEL -----------
[Time, P_data, dP_data, ddP_data, w2] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, 0, 1, [], qp_solver_type);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
    'color',[0 0.9 0], 'legend','$DMP^*_v$', 'plot3D',true, 'plot2D',true);

%% --------- Optimized DMP -> POS -----------
[Time, P_data, dP_data, ddP_data, w3] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, 1, 0, [], qp_solver_type);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
    'color',[0.72 0.27 1], 'legend','$DMP^*_p$', 'plot3D',true, 'plot2D',true);

%% NOTE about solver output
% "Solved inaccurate" is explained here https://osqp.org/docs/interfaces/status_values.html

%% ======== Plot Results ==========

%Plot demo in cartesian and unconstrained DMP result
poseDMP = ForwardKinematics(data{version}.Pos');
poseUnconstrainedDMP = ForwardKinematics(data{2}.Pos');

own_plots

%display whether limits are exceeded
min_joint_pos = (min(data{version}.Pos, [], 2) < actual_pos_lim(:,1))'
max_joint_pos = (max(data{version}.Pos, [], 2) > actual_pos_lim(:,2))'

min_joint_vel = (max(data{version}.Vel, [], 2) < actual_vel_lim(:,1))' 
max_joint_vel = (max(data{version}.Vel, [], 2) > actual_vel_lim(:,2))'

min_joint_acc = (max(data{version}.Accel, [], 2) < actual_accel_lim(:,1))'
max_joint_acc = (max(data{version}.Accel, [], 2) > actual_accel_lim(:,2))'



%% difference in joint velocity
%Calculate the difference between joint velocities in the unconstrained DMP and
%in the constrained DMP for the 1000 timesteps (1s total) with highest error. This way
%demonstrations with different lenght can be compared too.

d_jointVel = mean(maxk(abs(data{2}.Vel' - data{version}.Vel'), 1000), "all") %use max 1000 timesteps (total 1s time anywhere over the traj)

%% export
%Flip these joint signs so that output DMP runs on Franka2 in the lab
%wihtout sign flips, and flip signs back for Franka3 in the controller code
data{version}.Pos(1,:) = -data{version}.Pos(1,:);
data{version}.Pos(3,:) = -data{version}.Pos(3,:);
data{version}.Pos(5,:) = -data{version}.Pos(5,:);
data{version}.Pos(7,:) = -data{version}.Pos(7,:);
writematrix(data{version}.Pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat(bag,'_Opt_DMP_joint_',filename)))

data{version}.Vel(1,:) = -data{version}.Vel(1,:);
data{version}.Vel(3,:) = -data{version}.Vel(3,:);
data{version}.Vel(5,:) = -data{version}.Vel(5,:);
data{version}.Vel(7,:) = -data{version}.Vel(7,:);
writematrix(data{version}.Vel',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat(bag,'_Opt_DMP_joint_vel_',filename)))

%run forward kinematics again after sign flips so that cartesian trajectory
%can be plotted for comparison in control scripts
poseDMP = ForwardKinematics(data{version}.Pos');
writematrix(poseDMP,fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat(bag,'_Opt_DMP_pose_',filename)))
