%NOTE: modifications also to other files, e.g. in offlineGMPweigths number
%of DMP kernels are no longer reset to 30, number of points increased to
%400 from 200

clc;
close all;
clear;

addpath('utils/');
import_gmp_lib();
import_io_lib();

addpath('../../SupportScripts') %2 levels up due to Opt-DMP structure
addpath('../../Data/demos') %2 levels up due to Opt-DMP structure

%% Load training data
bag = 'A';
filename = '10l_bag_flip.csv';

%get demo_traj from this script
generateJointDemo

Timed = demo_traj.t;
Pd_data = demo_traj.pos;
dPd_data = demo_traj.vel;
ddPd_data = demo_traj.acc;

%% initialize and train GMP
train_method = 'LS';
N_kernels = 80; %Originally 30, difference in MSE but not noticable in plots
kernels_std_scaling = 4.0; %originally 1.5
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
p_max = [2.7437; 1.7628; 2.8973; -0.1518; 2.8065; 3.7525; 2.8973];
p_min = [-2.7437; -1.7628; -2.8973; -3.0421; -2.8065; 0.5445; -2.8973];
v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
a_max = [10; 7.5; 10; 10; 10; 10; 10];

pos_lim = [p_min p_max];
vel_lim = [-v_max v_max];
accel_lim = [-a_max a_max];

%use the actual robot limits for plotting
actual_pos_lim = pos_lim;
actual_vel_lim = vel_lim;
actual_accel_lim = accel_lim;

%use stricter limits than the actual ones for trajectory generation to have
%some margin
pos_lim = 0.98 * pos_lim;
vel_lim = 0.98 * vel_lim;
accel_lim = 0.98 * accel_lim;
%% ======== Generate trajectories ==========

qp_solver_type = 1; % matlab-quadprog:0 , osqp:1, Goldfarb-Idnani: 2

%% --------- Demo -----------
demo = struct('t',Timed, 'pos',Pd_data, 'vel',dPd_data, 'acc',ddPd_data, 'version','demo');

%% --------- Unconstrained DMP -----------
gmp2 = gmp.deepCopy();
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp2, tau, y0, yg);
unconstrained_DMP = struct('t',Time, 'pos',P_data, 'vel',dP_data, 'acc',ddP_data, 'version','UC-DMP');

%output unconstrained DMP for parameter comparison
writematrix(unconstrained_DMP.pos',fullfile('../../Data/trajectories',strcat(bag,'_UC_Opt_joint_',filename)))
writematrix(unconstrained_DMP.vel',fullfile('../../Data/trajectories',strcat(bag,'_UC_Opt_joint_vel_',filename)))

%% --------- Optimized DMP -> VEL -----------
%[Time, P_data, dP_data, ddP_data, w2] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, 0, 1, [], qp_solver_type);
%data{length(data)+1} = ...
%    struct('t',Time, 'pos',P_data, 'vel',dP_data, 'acc',ddP_data, 'version','$Opt-DMP$');

%% --------- Optimized DMP -> POS -----------
[Time, P_data, dP_data, ddP_data, w3] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, 1, 0, [], qp_solver_type);
res = struct('t',Time, 'pos',P_data, 'vel',dP_data, 'acc',ddP_data, 'version','Opt-DMP');

%% NOTE about solver output
% "Solved inaccurate" is explained here https://osqp.org/docs/interfaces/status_values.html

%% ======== Plot Results ==========

%Plot demo in cartesian and unconstrained DMP result
poseDMP = ForwardKinematics(res.pos');
poseUnconstrainedDMP = ForwardKinematics(unconstrained_DMP.pos');

%own_plots %OLD Opt-DMP specific script
generate_plots

%display whether limits are exceeded
min_joint_pos = (min(res.pos, [], 2) < actual_pos_lim(:,1))'
max_joint_pos = (max(res.pos, [], 2) > actual_pos_lim(:,2))'

min_joint_vel = (max(res.vel, [], 2) < actual_vel_lim(:,1))' 
max_joint_vel = (max(res.vel, [], 2) > actual_vel_lim(:,2))'

min_joint_acc = (max(res.acc, [], 2) < actual_accel_lim(:,1))'
max_joint_acc = (max(res.acc, [], 2) > actual_accel_lim(:,2))'


%% difference in joint velocity
%Calculate the difference between joint velocities in the unconstrained DMP and
%in the constrained DMP for the 1000 timesteps (1s total) with highest error. This way
%demonstrations with different lenght can be compared too.

d_jointVel = mean(maxk(abs(unconstrained_DMP.vel' - res.vel'), 1000), "all") %use max 1000 timesteps (total 1s time anywhere over the traj)


res_unflipped = res %for exporting manually for joint vel comparison plots

%% export
%Flip these joint signs so that output DMP runs on Franka2 in the lab
%wihtout sign flips, and flip signs back for Franka3 in the controller code
res.pos(1,:) = -res.pos(1,:);
res.pos(3,:) = -res.pos(3,:);
res.pos(5,:) = -res.pos(5,:);
res.pos(7,:) = -res.pos(7,:);
writematrix(res.pos',fullfile('../../Data/trajectories',strcat(bag,'_Opt_DMP_joint_',filename)))

res.vel(1,:) = -res.vel(1,:);
res.vel(3,:) = -res.vel(3,:);
res.vel(5,:) = -res.vel(5,:);
res.vel(7,:) = -res.vel(7,:);
writematrix(res.vel',fullfile('../../Data/trajectories',strcat(bag,'_Opt_DMP_joint_vel_',filename)))

%run forward kinematics again after sign flips so that cartesian trajectory
%can be plotted for comparison in control scripts
poseDMP = ForwardKinematics(res.pos');
writematrix(poseDMP,fullfile('../../Data/trajectories',strcat(bag,'_Opt_DMP_pose_',filename)))
