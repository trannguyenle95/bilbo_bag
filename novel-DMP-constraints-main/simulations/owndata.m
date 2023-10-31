clc;
close all;
clear;

addpath('utils/');
import_gmp_lib();
import_io_lib();

addpath('/home/erichannus/catkin_ws/src/SupportScripts/')
addpath('/home/erichannus/catkin_ws/src/Data/demos/')

%% Load training data
filename = '10l_bag_flip.csv';

version = 4 % 3 = vel, 4 = pos, select which optimization version to export

D = preprocess(filename, false, 0.00, 0.00, 0.00, 1, 'ori1', 0.38);
Dsmooth = smoothdata(D, 1, "gaussian",35); %smooth demo before calculating IK
Dsmooth(:,4:7) = Dsmooth(:,4:7) ./ sqrt(sum(Dsmooth(:,4:7).^2,2)); %Make sure quaternion is still unit
[q, jacobians] = InverseKinematics(Dsmooth);
demo_traj = generateDemo(q', 1/120);

Timed = demo_traj.t;
Pd_data = demo_traj.pos;
dPd_data = demo_traj.vel;
ddPd_data = demo_traj.acc;

%% initialize and train GMP
train_method = 'LS';
N_kernels = 200; %Originally 30, difference in MSE but not noticable in plots
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
pos_lim = 0.98 * pos_lim;
vel_lim = 0.98 * vel_lim;
accel_lim = 0.98 * accel_lim;

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

%save standard DMP traj
writematrix(data{2}.Pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('nominalDMP3_joint_',filename)))
writematrix(data{2}.Vel',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('nominalDMP3_joint_vel_',filename)))

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
label_font = 17;
ax_fontsize = 14;

ind = [1 2 3 4 5 6 7]; % choose DoFs to plot
for k=1:length(ind)
    i = ind(k);
    fig = figure;
    sgtitle(strcat('joint',int2str(i)))

    % plot joint positions
    ax = subplot(3,1,1);
    hold on;
    legend_ = {};
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Pos(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
        legend_ = [legend_ data{k}.legend];
    end
    axis tight;
    plot(ax.XLim, [actual_pos_lim(i,1) actual_pos_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [actual_pos_lim(i,2) actual_pos_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ylabel('pos [$rad$]', 'interpreter','latex', 'fontsize',label_font);
    legend(legend_, 'interpreter','latex', 'fontsize',17, 'Position',[0.2330 0.9345 0.5520 0.0294], 'Orientation', 'horizontal');
    ax.FontSize = ax_fontsize;
    hold off;

    % plot joint velocities
    ax = subplot(3,1,2);
    hold on;
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Vel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    plot(ax.XLim, [actual_vel_lim(i,1) actual_vel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [actual_vel_lim(i,2) actual_vel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ylabel('vel [$rad/s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

    % plot joint velocities
    ax = subplot(3,1,3);
    hold on;
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Accel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    plot(ax.XLim, [actual_accel_lim(i,1) actual_accel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [actual_accel_lim(i,2) actual_accel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ylabel('accel [$rad/s^2$]', 'interpreter','latex', 'fontsize',label_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;
end

% ======================================================

%Plot demo in cartesian and unconstrained DMP result
poseDMP = ForwardKinematics(data{2}.Pos'); %index2 = DMP without constraints
figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

plot(data{2}.Time,poseDMP(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_x')
plot(data{2}.Time,poseDMP(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_y')
plot(data{2}.Time,poseDMP(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(data{2}.Time,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_{qx}')
plot(data{2}.Time,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_{qy}')
plot(data{2}.Time,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_{qz}')
plot(data{2}.Time,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMP_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off


%Plot demo in cartesian and constrained DMP result
poseDMP = ForwardKinematics(data{3}.Pos'); %index 3 = vel optimization
figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

plot(data{3}.Time,poseDMP(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMPv_x')
plot(data{3}.Time,poseDMP(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMPv_y')
plot(data{3}.Time,poseDMP(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMPv_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(data{3}.Time,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMPv_{qx}')
plot(data{3}.Time,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMPv_{qy}')
plot(data{3}.Time,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMPv_{qz}')
plot(data{3}.Time,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMPv_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

%Plot demo in cartesian and constrained DMP result
poseDMP = ForwardKinematics(data{4}.Pos'); %index 4 = pos optimization
figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

plot(data{4}.Time,poseDMP(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMPp_x')
plot(data{4}.Time,poseDMP(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMPp_y')
plot(data{4}.Time,poseDMP(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMPp_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(data{4}.Time,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMPp_{qx}')
plot(data{4}.Time,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMPp_{qy}')
plot(data{4}.Time,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMPp_{qz}')
plot(data{4}.Time,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMPp_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off


%display whether limits are exceeded
min_joint_pos = (min(data{version}.Pos, [], 2) < actual_pos_lim(:,1))'
max_joint_pos = (max(data{version}.Pos, [], 2) > actual_pos_lim(:,2))'

min_joint_vel = (max(data{version}.Vel, [], 2) < actual_vel_lim(:,1))' 
max_joint_vel = (max(data{version}.Vel, [], 2) > actual_vel_lim(:,2))'

min_joint_acc = (max(data{version}.Accel, [], 2) < actual_accel_lim(:,1))'
max_joint_acc = (max(data{version}.Accel, [], 2) > actual_accel_lim(:,2))'



%% difference in joint velocity
%Calculate the difference between joint velocities in the demonstration and
%in the constrained DMP for the 1000 timesteps (1s total) with highest error. This way
%demonstrations with different lenght can be compared too.
demoVel = interpolate(data{1}.Vel', data{version}.Vel', false);
d_jointVel = mean(maxk(abs(demoVel - data{version}.Vel'), 1000), "all") %use max 1000 timesteps (total 1s time anywhere over the traj) 

%% export
%Flip these joint signs so that output DMP runs on Franka2 in the lab
%wihtout sign flips, and flip signs back for Franka3 in the controller code
data{version}.Pos(1,:) = -data{version}.Pos(1,:);
data{version}.Pos(3,:) = -data{version}.Pos(3,:);
data{version}.Pos(5,:) = -data{version}.Pos(5,:);
data{version}.Pos(7,:) = -data{version}.Pos(7,:);
writematrix(data{version}.Pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('DMP3_joint_',filename)))

data{version}.Vel(1,:) = -data{version}.Vel(1,:);
data{version}.Vel(3,:) = -data{version}.Vel(3,:);
data{version}.Vel(5,:) = -data{version}.Vel(5,:);
data{version}.Vel(7,:) = -data{version}.Vel(7,:);
writematrix(data{version}.Vel',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('DMP3_joint_vel_',filename)))

%run forward kinematics again after sign flips so that cartesian trajectory
%can be plotted for comparison in control scripts
poseDMP = ForwardKinematics(data{version}.Pos');
writematrix(poseDMP,fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('DMP3_pose_',filename)))
