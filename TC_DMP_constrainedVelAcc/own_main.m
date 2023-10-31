clearvars
clc
close all 

addpath('/home/erichannus/catkin_ws/src/SupportScripts/')
addpath('/home/erichannus/catkin_ws/src/Data/demos/')


% Create demonstration trajectory
filename = '10l_bag_flip.csv';
D = preprocess(filename, false, 0.00, 0.00, 0.00, 1, 'ori1', 0.38);
Dsmooth = smoothdata(D, 1, "gaussian",35); %smooth demo before calculating IK
Dsmooth(:,4:7) = Dsmooth(:,4:7) ./ sqrt(sum(Dsmooth(:,4:7).^2,2)); %Make sure quaternion is still unit
[q, jacobians] = InverseKinematics(Dsmooth);

%generate demo struct
demo_traj = generateDemo(q', 1/120);

% Nominal trajectory functions
dmp_params.D = 20;
dmp_params.K = dmp_params.D^2/4;
dmp_params.n_kernel = 800; %default 100
%Originally 100, better fit with 200 - with 100 some peaks from demo are reduced by smoothing
%BUT 200 can overfit and give acceleration that exceeds limits slightly,
%not enough smooothing effect in that case...
dmp_params.alpha_s = 1;
dmp = DMP(dmp_params,demo_traj);

% Simulation setup
sim_params.dt = 1/1000; %frequency of Franka robots
sim_params.T_max = 100; %maximum duration after slowdown due to constraints - original value 100s is much longer than demonstrated trajectory lengths so this value is fine

%Unconstrained DMP
nominalTraj = dmp.rollout(sim_params.dt);
writematrix(nominalTraj.pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('nominalDMP2_joint_',filename)))
writematrix(nominalTraj.vel',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('nominalDMP2_joint_vel_',filename)))


%limits are chosen to be worst case between Franka Panda and Franka
%Research 3
p_max = [2.7437; 1.7628; 2.8973; -0.1518; 2.8065; 3.7525; 2.8973]; %not used in sim, added for plotting
p_min = [-2.7437; -1.7628; -2.8973; -3.0421; -2.8065; 0.5445; -2.8973]; %not used in sim, added for plotting
sim_params.v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
sim_params.a_max = [10; 7.5; 10; 10; 10; 10; 10];
%note: vel and acc limits are scaled by 0.98 when passed to own_run_simulation.m for
%additional margin

% Scaling parameters
tc_params.nominal_tau = dmp.nominal_tau;
tc_params.eps = 1e-3;
tc_params.gamma_nominal = 1;
tc_params.gamma_a = 0.5; %originally 0.5, higher should reducde overshoots but increase runtime
tc_params.a_max = sim_params.a_max * 0.98; %added scaling here to have some margin as acc limits are not guaranteed
tc_params.v_max = sim_params.v_max * 0.98; %use same scaling for vel limits as for acc

% Create temporal coupling object
tc = TemporalCoupling(tc_params);

% Simulate
res = own_run_simulation(dmp,tc,sim_params);   

%print whether limits are exceeded
min_joint_pos = (min(res.ref_pos, [], 2) < p_min)'
max_joint_pos = (max(res.ref_pos, [], 2) > p_max)'

min_joint_vel = (min(res.ref_vel, [], 2) < -sim_params.v_max)' 
max_joint_vel = (max(res.ref_vel, [], 2) > sim_params.v_max)'

min_joint_acc = (min(res.ref_acc, [], 2) < -sim_params.a_max)'
max_joint_acc = (max(res.ref_acc, [], 2) > sim_params.a_max)'

slowdown = res.t(end) - demo_traj.t(end) %how much longer runtime (in seconds) does the constrained DMP have compared to the demonstrated trajectory

%calculate the corresponding cartesian trajectory for plotting
poseDMP = ForwardKinematics(res.ref_pos');

% Plot result
own_plot_result

%Flip these joint signs so that output DMP runs on Franka2 in the lab
%wihtout sign flips, and flip signs back for Franka3 in the controller code
res.ref_pos(1,:) = -res.ref_pos(1,:);
res.ref_pos(3,:) = -res.ref_pos(3,:);
res.ref_pos(5,:) = -res.ref_pos(5,:);
res.ref_pos(7,:) = -res.ref_pos(7,:);
writematrix(res.ref_pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('DMP2_joint_',filename)))

res.ref_vel(1,:) = -res.ref_vel(1,:);
res.ref_vel(3,:) = -res.ref_vel(3,:);
res.ref_vel(5,:) = -res.ref_vel(5,:);
res.ref_vel(7,:) = -res.ref_vel(7,:);
writematrix(res.ref_vel',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('DMP2_joint_vel_',filename)))

%run forward kinematics again after sign flips so that cartesian trajectory
%can be plotted for comparison in control scripts
poseDMP = ForwardKinematics(res.ref_pos');
writematrix(poseDMP,fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('DMP2_pose_',filename)))
