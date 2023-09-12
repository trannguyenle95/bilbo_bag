clearvars
clc
close all 

% Create demonstration trajectory
filename = 'sack_from_bag2.csv';
D = preprocess(filename, false, 2);

q = InverseKinematics2(D)';

poseIKFK = ForwardKinematics2(q');

demo_traj = generateDemo(q, 1/120);


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
sim_params.T_max = 100;

nominalTraj = dmp.rollout(sim_params.dt);


%WORST CASE COMMON LIMITS
p_max = [2.7437; 1.7628; 2.8973; -0.1518; 2.8065; 3.7525; 2.8973]; %not used in sim, added for PLOTTING
p_min = [-2.7437; -1.7628; -2.8973; -3.0421; -2.8065; 0.5445; -2.8973]; %not used in sim, added for PLOTTING
sim_params.v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
sim_params.a_max = [10; 7.5; 10; 10; 10; 10; 10];

% Scaling parameters
tc_params.nominal_tau = dmp.nominal_tau;
tc_params.eps = 1e-3;
tc_params.gamma_nominal = 1;
tc_params.gamma_a = 3; %originally 0.5, higher should reducde overshoots but increase runtime
tc_params.a_max = sim_params.a_max; 
tc_params.v_max = sim_params.v_max;

% Create temporal coupling objects
tc = {};
tc{1} = TemporalCoupling(tc_params);
%tc{2} = NoTemporalCoupling(tc_params);
%tc{3} = TemporalCouplingNoPotentials(tc_params);

% Simulate
nSim = 5;
res = cell(length(tc),1);

pose = cell(length(tc),1);
for d = 1:nSim
    disp(['Running simulation ' num2str(d) '/' num2str(length(tc)) '...']);
    res{d} = run_simulation(dmp,tc{1},sim_params);  
    pose{d} = ForwardKinematics2(res{d}.ref_pos');
end

figure('Name','x')
hold on
plot(res{1}.t,pose{1}(:,1),'LineWidth',3)
plot(res{2}.t,pose{2}(:,1),'LineWidth',2.5)
plot(res{3}.t,pose{3}(:,1),'LineWidth',2)
plot(res{4}.t,pose{4}(:,1),'LineWidth',1.5)
plot(res{5}.t,pose{5}(:,1),'LineWidth',1)

hold off
