clearvars
clc
close all 

addpath('../SupportScripts')
addpath('../Data/demos')


% Create demonstration trajectory

bag = 'A';
filename = '10l_bag_flip.csv';
gamma_a = 0.65; %originally 0.5, higher should reducde overshoots but increase runtime

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

% Simulation setup
sim_params.dt = 1/1000; %frequency of Franka robots
sim_params.T_max = 100; %maximum duration after slowdown due to constraints - original value 100s is much longer than demonstrated trajectory lengths so this value is fine

%Unconstrained DMP
unconstrained_DMP = dmp.rollout(sim_params.dt, 1);

%output unconstrained DMP for parameter comparison
writematrix(unconstrained_DMP.pos',fullfile('../Data/trajectories',strcat(bag,'_UC_TC_DMP_joint_',filename)))
writematrix(unconstrained_DMP.vel',fullfile('../Data/trajectories',strcat(bag,'_UC_TC_DMP_joint_vel_',filename)))

%limits are chosen to be worst case between Franka Panda and Franka
%Research 3
p_max = [2.7437; 1.7628; 2.8973; -0.1518; 2.8065; 3.7525; 2.8973]; %not used in sim, added for plotting
p_min = [-2.7437; -1.7628; -2.8973; -3.0421; -2.8065; 0.5445; -2.8973]; %not used in sim, added for plotting
v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
a_max = [10; 7.5; 10; 10; 10; 10; 10];
sim_params.v_max = v_max;
sim_params.a_max = a_max;


%note: vel and acc limits are scaled by 0.98 when passed to own_run_simulation.m for
%additional margin

% Scaling parameters
tc_params.nominal_tau = dmp.nominal_tau;
tc_params.eps = 1e-3;
tc_params.gamma_nominal = 1;
tc_params.gamma_a = gamma_a; %originally 0.5, higher should reducde overshoots but increase runtime
tc_params.a_max = sim_params.a_max * 0.98; %added scaling here to have some margin as acc limits are not guaranteed
tc_params.v_max = sim_params.v_max * 0.98; %use same scaling for vel limits as for acc

% Create temporal coupling object
tc = TemporalCoupling(tc_params);

% Simulate
res = own_run_simulation(dmp,tc,sim_params);   

%print whether limits are exceeded
min_joint_pos = (min(res.pos, [], 2) < p_min)'
max_joint_pos = (max(res.pos, [], 2) > p_max)'

min_joint_vel = (min(res.vel, [], 2) < -sim_params.v_max)' 
max_joint_vel = (max(res.vel, [], 2) > sim_params.v_max)'

min_joint_acc = (min(res.acc, [], 2) < -sim_params.a_max)'
max_joint_acc = (max(res.acc, [], 2) > sim_params.a_max)'

slowdown = (res.t(end) - demo_traj.t(end)) / demo_traj.t(end)*100; %how much longer runtime (%) does the constrained DMP have compared to the demonstrated trajectory
disp(strcat('slowdown: ', num2str(slowdown), ' %'))

%calculate the corresponding cartesian trajectory for plotting
poseUnconstrainedDMP = ForwardKinematics(unconstrained_DMP.pos');
poseDMP = ForwardKinematics(res.pos');

% Plot result
res.version = "TC-DMP"
generate_plots

%%%
% Tau vs time (high tau means that slowed down motion is executed durnig these timesteps)
figure('Name','Tau')
plot(res.t,res.tau,'b','LineWidth',2.5), hold on
set(gca,'FontSize',16)
ylim([0.9*dmp.nominal_tau 1.05*max(res.tau)])
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
ylabel('$\tau$', 'Interpreter','latex','Fontsize',30)



res_unflipped = res

%Flip these joint signs so that output DMP runs on Franka2 in the lab
%wihtout sign flips, and flip signs back for Franka3 in the controller code
res.pos(1,:) = -res.pos(1,:);
res.pos(3,:) = -res.pos(3,:);
res.pos(5,:) = -res.pos(5,:);
res.pos(7,:) = -res.pos(7,:);
writematrix(res.pos',fullfile('../Data/trajectories',strcat(bag,'_TC_DMP_joint_',filename)))

res.vel(1,:) = -res.vel(1,:);
res.vel(3,:) = -res.vel(3,:);
res.vel(5,:) = -res.vel(5,:);
res.vel(7,:) = -res.vel(7,:);
writematrix(res.vel',fullfile('../Data/trajectories',strcat(bag,'_TC_DMP_joint_vel_',filename)))

res.acc(1,:) = -res.acc(1,:);
res.acc(3,:) = -res.acc(3,:);
res.acc(5,:) = -res.acc(5,:);
res.acc(7,:) = -res.acc(7,:);
writematrix(res.acc',fullfile('../Data/trajectories',strcat(bag,'_TC_DMP_joint_acc_',filename)))

%run forward kinematics again after sign flips so that cartesian trajectory
%can be plotted for comparison in control scripts
poseDMP = ForwardKinematics(res.pos');
writematrix(poseDMP,fullfile('../Data/trajectories',strcat(bag,'_TC_DMP_pose_',filename)))
