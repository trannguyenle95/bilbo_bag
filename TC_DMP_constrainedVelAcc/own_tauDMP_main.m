clearvars
clc
close all 

addpath('../SupportScripts')
addpath('../Data/demos')

% Create demonstration trajectory
bag = 'A';
filename = '10l_bag_flip.csv';
tau_scale = 3.9;

%get demo_traj from this script
generateJointDemo

% Nominal trajectory functions
dmp_params.D = 20;
dmp_params.K = dmp_params.D^2/4;
dmp_params.n_kernel = 800; %default 100
%Originally 100, better fit with 200 - with 100 some peaks from demo are reduced by smoothing
%BUT 200 can overfit and give acceleration that exceeds limits slightly,
%not enough smooothing effect in that case...
dmp_params.alpha_s = 1;
dmp = DMP(dmp_params,demo_traj);


dt = 1/1000; %frequency of Franka robots

%Unconstrained DMP
unconstrained_DMP = dmp.rollout(dt, 1);

%DMP constrained by manually setting tau
res = dmp.rollout(dt, tau_scale);

%limits for plotting and checking if they are exceeded
p_max = [2.7437; 1.7628; 2.8973; -0.1518; 2.8065; 3.7525; 2.8973]; %not used in sim, added for plotting
p_min = [-2.7437; -1.7628; -2.8973; -3.0421; -2.8065; 0.5445; -2.8973]; %not used in sim, added for plotting
v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
a_max = [10; 7.5; 10; 10; 10; 10; 10];

%print whether limits are exceeded
min_joint_pos = (min(res.pos, [], 2) < p_min)'
max_joint_pos = (max(res.pos, [], 2) > p_max)'

min_joint_vel = (min(res.vel, [], 2) < -v_max)' 
max_joint_vel = (max(res.vel, [], 2) > v_max)'

min_joint_acc = (min(res.acc, [], 2) < -a_max)'
max_joint_acc = (max(res.acc, [], 2) > a_max)'

slowdown = (res.t(end) - demo_traj.t(end)) / demo_traj.t(end)*100; %how much longer runtime (%) does the constrained DMP have compared to the demonstrated trajectory
disp(strcat('slowdown: ', num2str(slowdown), ' %'))


%calculate the corresponding cartesian trajectory for plotting
poseUnconstrainedDMP = ForwardKinematics(unconstrained_DMP.pos');
poseDMP = ForwardKinematics(res.pos');

% Plot result
res.version = "tau-DMP"
generate_plots

res_unflipped = res

%Flip these joint signs so that output DMP runs on Franka2 in the lab
%wihtout sign flips, and flip signs back for Franka3 in the controller code
res.pos(1,:) = -res.pos(1,:);
res.pos(3,:) = -res.pos(3,:);
res.pos(5,:) = -res.pos(5,:);
res.pos(7,:) = -res.pos(7,:);
writematrix(res.pos',fullfile('../Data/trajectories',strcat(bag,'_tau_DMP_joint_',filename)))

res.vel(1,:) = -res.vel(1,:);
res.vel(3,:) = -res.vel(3,:);
res.vel(5,:) = -res.vel(5,:);
res.vel(7,:) = -res.vel(7,:);
writematrix(res.vel',fullfile('../Data/trajectories',strcat(bag,'_tau_DMP_joint_vel_',filename)))

res.acc(1,:) = -res.acc(1,:);
res.acc(3,:) = -res.acc(3,:);
res.acc(5,:) = -res.acc(5,:);
res.acc(7,:) = -res.acc(7,:);
writematrix(res.acc',fullfile('../Data/trajectories',strcat(bag,'_tau_DMP_joint_acc_',filename)))

%run forward kinematics again after sign flips so that cartesian trajectory
%can be plotted for comparison in control scripts
poseDMP = ForwardKinematics(res.pos');
writematrix(poseDMP,fullfile('../Data/trajectories',strcat(bag,'_tau_DMP_pose_',filename)))
