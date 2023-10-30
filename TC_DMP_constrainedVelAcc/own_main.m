clearvars
clc
close all 

addpath('/home/erichannus/catkin_ws/src/MatlabPreProcess/')
addpath('/home/erichannus/catkin_ws/src/Data/demos/')


% Create demonstration trajectory
%filename = 'sack_from_bag2.csv';
%D = preprocess(filename, false, 0.57, 0, -0.20, 1); %other processing of x distance to avoid collisions... >>  if larger range of x allowed then robots need to be moved closer potentially

filename = '10l_bag_flip.csv';
D = preprocess(filename, false, 0.00, 0.00, 0.00, 1, 'ori2', 0.38);
%D = preprocess(filename, false, 0.05, -0.05, 0.05, 1);

Dsmooth = smoothdata(D, 1, "gaussian",35); %still noisy after IK
Dsmooth(:,4:7) = Dsmooth(:,4:7) ./ sqrt(sum(Dsmooth(:,4:7).^2,2)); %Make sure quaternion is still unit

[q, jacobians] = InverseKinematics(Dsmooth);

poseIKFK = ForwardKinematics(q);

%q = smoothdata(q, 2, "gaussian",10);



%generate demo struct
delta_t = 1/120;
demo_traj.t = delta_t * [0:1:length(q')-1];
demo_traj.pos = q';
demo_traj.vel = [diff(q',1,2)/delta_t zeros(size(demo_traj.pos,1), 1)];
demo_traj.acc = [diff(demo_traj.vel,1,2)/delta_t zeros(size(demo_traj.pos,1), 1)];


% Nominal trajectory functions
dmp_params.D = 20;
dmp_params.K = dmp_params.D^2/4;
dmp_params.n_kernel = 200; %default 100
%Originally 100, better fit with 200 - with 100 some peaks from demo are reduced by smoothing
%BUT 200 can overfit and give acceleration that exceeds limits slightly,
%not enough smooothing effect in that case...
dmp_params.alpha_s = 1;
dmp = DMP(dmp_params,demo_traj);

% Simulation setup
sim_params.dt = 1/1000; %NOTE: this is actually correct rate for Frankas too!
%sim_params.dt = 1/120; %Test with 1/120 to fit DMP with OptiTrack frequency
sim_params.T_max = 100;
%sim_params.kv = 10;
%sim_params.kp = sim_params.kv^2/4;

nominalTraj = dmp.rollout(sim_params.dt); %ROLLOUT only for nominal traj? Then use DMP.step in simulation script to get constrained DMP result?
%sim_params.a_max = 0.5*[max(abs(nominalTraj.acc(1,:))); max(abs(nominalTraj.acc(2,:)))];
%sim_params.v_max = inf*[1 1]'; %NOTE: change this if I want to have v limit from the start
%sim_params.v_max = [0.15 0.15]';
%sim_params.v_max = 0.5*[max(abs(nominalTraj.vel(1,:))); max(abs(nominalTraj.vel(2,:)))];
%sim_params.s_v_max_online = 0.65;
%sim_params.v_max_online = [inf;0.15]; %leave defined here for plotting, but disable in sim
% sim_params.v_max_online = inf*[1 1]'; %Commented away in original code,
% can be uncommented to remove limit


%WORST CASE COMMON LIMITS
p_max = [2.7437; 1.7628; 2.8973; -0.1518; 2.8065; 3.7525; 2.8973]; %not used in sim, added for PLOTTING
p_min = [-2.7437; -1.7628; -2.8973; -3.0421; -2.8065; 0.5445; -2.8973]; %not used in sim, added for PLOTTING
sim_params.v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
sim_params.a_max = [10; 7.5; 10; 10; 10; 10; 10];
%NOTE: scale vel and acc limits by 0.98 in run_simulation.m

% Scaling parameters
tc_params.nominal_tau = dmp.nominal_tau;
tc_params.eps = 1e-3;
tc_params.gamma_nominal = 1;
tc_params.gamma_a = 0.5; %originally 0.5, higher should reducde overshoots but increase runtime
tc_params.a_max = sim_params.a_max * 0.98; %NOTE added scaling here to have some margin as vel/acc limits are not guaranteed
tc_params.v_max = sim_params.v_max * 0.98; %NOTE added scaling here to have some margin as vel/acc limits are not guaranteed
%alternativley gamma_a can be increased to get rid of small overshoots, but
%trajectory slowdown experimentally seems to be lower when adding small
%margin by scaling 0.98

% Create temporal coupling objects
tc = {};
tc{1} = TemporalCoupling(tc_params);
%tc{2} = NoTemporalCoupling(tc_params);
%tc{3} = TemporalCouplingNoPotentials(tc_params);

%TODO: clean up this + plotting, no need to use cells when only tc{1} is
%used
% Simulate
nSim = length(tc);
res = cell(length(tc),1);
for d = 1:nSim
    disp(['Running simulation ' num2str(d) '/' num2str(length(tc)) '...']);
    res{d} = own_run_simulation(dmp,tc{d},sim_params);   
end

%plot whether limits are exceeded
min_joint_pos = (min(res{1}.ref_pos, [], 2) < p_min)'
max_joint_pos = (max(res{1}.ref_pos, [], 2) > p_max)'

min_joint_vel = (min(res{1}.ref_vel, [], 2) < -sim_params.v_max)' 
max_joint_vel = (max(res{1}.ref_vel, [], 2) > sim_params.v_max)'

min_joint_acc = (min(res{1}.ref_acc, [], 2) < -sim_params.a_max)'
max_joint_acc = (max(res{1}.ref_acc, [], 2) > sim_params.a_max)'

slowdown = res{1}.t(end) - demo_traj.t(end)

%first do this for plotting
poseDMP = ForwardKinematics(res{1}.ref_pos');

% Plot result
own_plot_result

%also do this for both robots and let controller do sign flips
%NOTE: later change controller to use Franka3 as default and these values
%would not need to be flipped? Note jonint7 is not just a flip but also
%some offset (see controller code)
res{1}.ref_pos(1,:) = -res{1}.ref_pos(1,:);
res{1}.ref_pos(3,:) = -res{1}.ref_pos(3,:);
res{1}.ref_pos(5,:) = -res{1}.ref_pos(5,:);
res{1}.ref_pos(7,:) = -res{1}.ref_pos(7,:);

%Run Forward Kinematics for plotting and playback with cartesian control
%previously and do this again so that pose with flipped joints is saved for pose
%trajectory
poseDMP = ForwardKinematics(res{1}.ref_pos');

%save joint angles to file
%writematrix(res{1}.ref_pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('joint_',filename)))
%writematrix(poseDMP,fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('pose_',filename)))

%use general filename for to easily switch between demos when generating
%from new files above
writematrix(res{1}.ref_pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories','joint_demoDMP.csv'))
writematrix(poseDMP,fullfile('/home/erichannus/catkin_ws/src/Data/trajectories','pose_demoDMP.csv'))

res{1}.ref_vel(1,:) = -res{1}.ref_vel(1,:);
res{1}.ref_vel(3,:) = -res{1}.ref_vel(3,:);
res{1}.ref_vel(5,:) = -res{1}.ref_vel(5,:);
res{1}.ref_vel(7,:) = -res{1}.ref_vel(7,:);

%writematrix(res{1}.ref_vel',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('joint_vel_',filename)))
writematrix(res{1}.ref_vel',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories','joint_vel_demoDMP.csv'))