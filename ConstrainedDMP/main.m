clearvars
clc
close all 

% Create demonstration trajectory
%load('trajectories/g')

%D = readmatrix(strcat('demos/','BagFlip_q.csv'))';


%filename = 'BagFlip.csv';
%D = preprocess_1hand(filename, false);
filename = 'BagFlip.csv';
D = preprocess(filename, false, 2);

%D = smoothdata(D, 1, "gaussian",50); %still noisy after IK

q = InverseKinematics(D)';

%q = smoothdata(q, 2, "gaussian",10);

demo_traj = generateDemo(q, 1/120);
%demo_traj = generateDemo('sack2_pre_pyth_ik_mat_q.csv', 1/120)

%demo_traj = generateDemo('sack2_pre_pyth_ik_mat_q.csv', 1/120)

%NOTE: this dmeo too large matrix to handle??
% TODO: how to replay on robots? Need to go from tau etc back to setpoints
% for each timestep 1 ms...

% load('trajectories/omega')
% load('trajectories/z')
% load('trajectories/s')
% load('trajectories/psi')

% Nominal trajectory functions
dmp_params.D = 20;
dmp_params.K = dmp_params.D^2/4;
dmp_params.n_kernel = 100; %default 100
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

%FOR FRANKA PANDA (TODO RESEARCH 3 + toggle with input argument)
p_max = [2.8973; 1.7628; 2.8973; -0.0698; 2.8973; 3.7525; 2.8973]; %not used in sim, added for PLOTTING
p_min = [-2.8973; -1.7628; -2.8973; -3.0718; -2.8973; -0.0175; -2.8973]; %not used in sim, added for PLOTTING
sim_params.v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
sim_params.a_max = [15; 7.5; 10; 12.5; 15; 20; 20];

% Scaling parameters
tc_params.nominal_tau = dmp.nominal_tau;
tc_params.eps = 1e-3;
tc_params.gamma_nominal = 1;
tc_params.gamma_a = 1; %originally 0.5, higher should reducde overshoots but increase runtime
tc_params.a_max = sim_params.a_max; 
tc_params.v_max = sim_params.v_max;

% Create temporal coupling objects
tc = {};
tc{1} = TemporalCoupling(tc_params);
%tc{2} = NoTemporalCoupling(tc_params);
%tc{3} = TemporalCouplingNoPotentials(tc_params);

% Simulate
nSim = length(tc);
res = cell(length(tc),1);
for d = 1:nSim
    disp(['Running simulation ' num2str(d) '/' num2str(length(tc)) '...']);
    res{d} = run_simulation(dmp,tc{d},sim_params);   
end


min_joint_pos = (min(res{1}.ref_pos, [], 2) < p_min)'
max_joint_pos = (max(res{1}.ref_pos, [], 2) > p_max)'

min_joint_vel = (min(res{1}.ref_vel, [], 2) < -sim_params.v_max)' 
max_joint_vel = (max(res{1}.ref_vel, [], 2) > sim_params.v_max)'

min_joint_acc = (min(res{1}.ref_acc, [], 2) < -sim_params.a_max)'
max_joint_acc = (max(res{1}.ref_acc, [], 2) > sim_params.a_max)'

%Run Forward Kinematics for plotting and playback with cartesian control
poseDMP = ForwardKinematics(res{1}.ref_pos');

%q_franka3 = res{1}.ref_pos';
%q_franka3(:,1) = -q_franka3(:,1);
%q_franka3(:,3) = -q_franka3(:,3);
%q_franka3(:,5) = -q_franka3(:,5);
%q_franka3(:,7) = -q_franka3(:,7);

% Plot result
plot_result

%save joint angles to file
writematrix(res{1}.ref_pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('joint_',filename)))
%writematrix(q_franka3,fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('joint_',filename)))
writematrix(poseDMP,fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('pose_',filename)))
