clc;
close all;
clear;

addpath('utils/');
import_gmp_lib();
import_io_lib();

addpath('own_code/') %ADDED

%% Load training data
% fid = FileIO('data/pos_data.bin', FileIO.in);
% Timed = fid.read('Timed');
% Pd_data = fid.read('Pd_data');
% dPd_data = fid.read('dPd_data');
% ddPd_data = fid.read('ddPd_data');
% fid.close();
%NOTE ^^ skip FileIO above and just import own data in shape of Pd, dPd and
%ddPd data!
%Each is an array of doubles, size DOF*timesteps

filename = '2h_flip2.csv';
D = preprocess(filename, false, 0.57, 0, 0, 1);
%D = preprocess(filename, false, 0.60, 0, 0, 1);
%D = preprocess(filename, false, 0.05, -0.05, 0.05, 1);
Dsmooth = smoothdata(D, 1, "gaussian",35); %still noisy after IK
Dsmooth(:,4:7) = Dsmooth(:,4:7) ./ sqrt(sum(Dsmooth(:,4:7).^2,2)); %Make sure quaternion is still unit
[q, jacobians] = InverseKinematics2(Dsmooth);
demo_traj = generateDemo(q', 1/120);

Timed = demo_traj.t;
Pd_data = demo_traj.pos;
dPd_data = demo_traj.vel;
ddPd_data = demo_traj.acc;

%Ts = Timed(2)-Timed(1); %NOTE: why isn't TS used?

%% initialize and train GMP
train_method = 'LS'; %NOTE: I can also try LWR
N_kernels = 100; %Originally 30, difference in MSE but not noticable in plots
kernels_std_scaling = 1.5;
n_dof = size(Pd_data,1);
gmp = GMP(n_dof, N_kernels, kernels_std_scaling); %NOTE: there is also GMPo for Cartesian pose
tic
offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
offline_train_mse
toc

gmp.setScaleMethod( TrajScale_Prop(n_dof) );

taud = Timed(end);
yd0 = gmp.getYd(0); %Pd_data(1);
ygd = gmp.getYd(1); %Pd_data(end);

%kt = 1.5; % temporal scaling %NOTE <<<<<<<<<<<<<<<<< CHECK THIS originally 1.5
kt = 1; %Disables temporal scaling of demo traj
%ks = diag([1 1 1]); % spatial scaling
tau = taud/kt;
y0 = yd0 + 0;
% yg = ks*(ygd - yd0) + y0;
% yg = ygd + [0.1; -0.1; 0.23]; view_ = [171.5301, -2.3630];
%yg = ygd + [0.7; -0.7; 0.05];  view_ = [171.9421, -3.0690]; %NOTE: SHIFT GOAL - original
yg = ygd; %Rollout goal is now same as demo goal!


%% ======== Limits ==========
%            lower limit     upper limit
%pos_lim = [[-1.2 -1.2 0.22]' [1.2 1.2 0.5]'];
%vel_lim = repmat([-0.3 0.3], 3,1);  % lower and upper limit, same for all DoFs
%accel_lim = repmat([-0.4 0.4], 3,1);

pos_lim = [[-2.7437; -1.7628; -2.8973; -3.0421; -2.8065; 0.5445; -2.8973] [2.7437; 1.7628; 2.8973; -0.1518; 2.8065; 3.7525; 2.8973]];
vel_lim = [-[2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100] [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100]];  % lower and upper limit, same for all DoFs
accel_lim = [-[10; 7.5; 10; 10; 10; 10; 10] [10; 7.5; 10; 10; 10; 10; 10]];

actual_pos_lim = pos_lim;
actual_vel_lim = vel_lim;
actual_accel_lim = accel_lim;

%use stricter limits than the actual ones
pos_lim = 0.98 * pos_lim;
vel_lim = 0.98 * vel_lim;
accel_lim = 0.98 * accel_lim;

%% ======== Generate trajectories ==========

qp_solver_type = 1; % matlab-quadprog:0 , osqp:1, Goldfarb-Idnani: 2

data = {};

%% --------- Demo -----------
data{length(data)+1} = ...
    struct('Time',Timed/kt, 'Pos',Pd_data, 'Vel',dPd_data*kt, 'Accel',ddPd_data*kt^2, 'linestyle','--', ...
    'color',[0.7 0.7 0.7], 'legend','demo', 'plot3D',true, 'plot2D',true);

%NOTE: Figure out if dPd and ddPd are actual vel/acc or just diff... <<<
%seems to be actual vel/acc same as in previous DMP codes

%% --------- Proportional scaling -----------
gmp2 = gmp.deepCopy();
gmp2.setScaleMethod( TrajScale_Prop(n_dof) );
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp2, tau, y0, yg);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
    'color','blue', 'legend','DMP', 'plot3D',true, 'plot2D',true);

%NOTE: this produces normal DMP (blue in paper)

%% --------- Optimized DMP -> VEL -----------
[Time, P_data, dP_data, ddP_data] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, 0, 1, [], qp_solver_type); %getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, false, true);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
    'color',[0 0.9 0], 'legend','$DMP^*_v$', 'plot3D',true, 'plot2D',true);

%% --------- Optimized DMP -> POS -----------
[Time, P_data, dP_data, ddP_data] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, 1, 0, [], qp_solver_type); %getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, true, false);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
    'color',[0.72 0.27 1], 'legend','$DMP^*_p$', 'plot3D',true, 'plot2D',true);

%% NOTE about solver output
% "Solved inaccurate" is explained here https://osqp.org/docs/interfaces/status_values.html

%% ======== Plot Results ==========

%title_ = {'x coordinate', 'y coordinate', 'z coordinate'};
label_font = 17;
ax_fontsize = 14;

ind = [1 2 3 4 5 6 7]; % choose DoFs to plot, e.g. [1 2 3] for [x y z]
for k=1:length(ind)
    i = ind(k);
    fig = figure;

    ax = subplot(3,1,1);
    hold on;
    % plot position trajectory
    legend_ = {};
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Pos(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
        legend_ = [legend_ data{k}.legend];
    end
    axis tight;
    % plot start and final positions - NOTE: commented out as it is not
    % needed at the moment
    %plot(0, y0(i), 'LineWidth', 4, 'LineStyle','none', 'Color','green','Marker','o', 'MarkerSize',10);
    %plot(tau, yg(i), 'LineWidth', 4, 'LineStyle','none', 'Color','red','Marker','x', 'MarkerSize',10);
    %plot(tau, ygd(i), 'LineWidth', 4, 'LineStyle','none', 'Color','magenta','Marker','x', 'MarkerSize',10);
    % plot bounds
    plot(ax.XLim, [actual_pos_lim(i,1) actual_pos_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [actual_pos_lim(i,2) actual_pos_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    % labels, title ...
    %ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',label_font);
    ylabel('pos [$rad$]', 'interpreter','latex', 'fontsize',label_font);
%     title(title_{i}, 'interpreter','latex', 'fontsize',18);
    legend(legend_, 'interpreter','latex', 'fontsize',17, 'Position',[0.2330 0.9345 0.5520 0.0294], 'Orientation', 'horizontal');
    ax.FontSize = ax_fontsize;
    hold off;

    ax = subplot(3,1,2);
    hold on;
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Vel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    % plot bounds
    plot(ax.XLim, [actual_vel_lim(i,1) actual_vel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [actual_vel_lim(i,2) actual_vel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    %ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',label_font);
    ylabel('vel [$rad/s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

    ax = subplot(3,1,3);
    hold on;
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Accel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    % plot bounds
    plot(ax.XLim, [actual_accel_lim(i,1) actual_accel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [actual_accel_lim(i,2) actual_accel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    %ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',label_font);
    ylabel('accel [$rad/s^2$]', 'interpreter','latex', 'fontsize',label_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

end

% ======================================================

%Plot demo in cartesian and constrained DMP result
poseDMP = ForwardKinematics2(data{2}.Pos'); %index2 = DMP without constraints
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
poseDMP = ForwardKinematics2(data{3}.Pos'); %index 3 = vel optimization
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
poseDMP = ForwardKinematics2(data{4}.Pos'); %index 4 = pos optimization
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




version = 3 % 3 = vel, 4 = pos;

%display whether limits are exceeded
min_joint_pos = (min(data{version}.Pos, [], 2) < actual_pos_lim(:,1))'
max_joint_pos = (max(data{version}.Pos, [], 2) > actual_pos_lim(:,2))'

min_joint_vel = (max(data{version}.Vel, [], 2) < actual_vel_lim(:,1))' 
max_joint_vel = (max(data{version}.Vel, [], 2) > actual_vel_lim(:,2))'

min_joint_acc = (max(data{version}.Accel, [], 2) < actual_accel_lim(:,1))'
max_joint_acc = (max(data{version}.Accel, [], 2) > actual_accel_lim(:,2))'

%also do this for both robots and let controller do sign flips
%NOTE: later change controller to use Franka3 as default and these values
%would not need to be flipped? Note jonint7 is not just a flip but also
%some offset (see controller code)
data{version}.Pos(1,:) = -data{version}.Pos(1,:);
data{version}.Pos(3,:) = -data{version}.Pos(3,:);
data{version}.Pos(5,:) = -data{version}.Pos(5,:);
data{version}.Pos(7,:) = -data{version}.Pos(7,:);

%Run Forward Kinematics for plotting and playback with cartesian control
%do this again so that pose with flipped joints is saved for pose
%trajectory
poseDMP = ForwardKinematics2(data{version}.Pos');

%save joint angles to file
writematrix(data{version}.Pos',fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('joint_',filename)))
writematrix(poseDMP,fullfile('/home/erichannus/catkin_ws/src/Data/trajectories',strcat('pose_',filename)))
%writematrix(data{3}.Pos',fullfile('output',strcat('joint_',filename)))
%writematrix(poseDMP,fullfile('output',strcat('pose_',filename)))