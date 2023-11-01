clearvars
clc
close all 

addpath('/home/erichannus/catkin_ws/src/Data/trajectories/')
TC_DMPpos = csvread('B_UC_TC_DMP_joint_10l_bag_flip.csv');
TC_DMPvel = csvread('B_UC_TC_DMP_joint_vel_10l_bag_flip.csv');
Opt_DMPpos = csvread('B_UC_Opt_joint_10l_bag_flip.csv');
Opt_DMPvel = csvread('B_UC_Opt_joint_vel_10l_bag_flip.csv');
%NOTE: change bag width below to match bag!

addpath('/home/erichannus/catkin_ws/src/SupportScripts/')
addpath('/home/erichannus/catkin_ws/src/Data/demos/')

% Create demonstration trajectory
filename = '10l_bag_flip.csv';
D = preprocess(filename, false, 0.00, 0.00, 0.00, 1, 'ori1', 0.37); %NOTE: CHANGE ACCORDING TO BAG
Dsmooth = smoothdata(D, 1, "gaussian",50); %smooth demo before calculating IK
Dsmooth(:,4:7) = Dsmooth(:,4:7) ./ sqrt(sum(Dsmooth(:,4:7).^2,2)); %Make sure quaternion is still unit
[q, jacobians] = InverseKinematics(Dsmooth);

% D = preprocess(filename, false, 0.00, 0.00, 0.00, 1, 'ori1', 0.38);
% Dsmooth = smoothdata(D, 1, "gaussian",35); %smooth demo before calculating IK
% Dsmooth(:,4:7) = Dsmooth(:,4:7) ./ sqrt(sum(Dsmooth(:,4:7).^2,2)); %Make sure quaternion is still unit
% [q, jacobians] = InverseKinematics(Dsmooth);
% 
% Dsmooth2 = smoothdata(D, 1, "gaussian",50); %smooth demo before calculating IK
% Dsmooth2(:,4:7) = Dsmooth2(:,4:7) ./ sqrt(sum(Dsmooth2(:,4:7).^2,2)); %Make sure quaternion is still unit
% 
% Dsmooth3 = smoothdata(D, 1, "gaussian",65); %smooth demo before calculating IK
% Dsmooth3(:,4:7) = Dsmooth3(:,4:7) ./ sqrt(sum(Dsmooth3(:,4:7).^2,2)); %Make sure quaternion is still unit

% for DOF = 1:7
%     figure('Name',strcat('DOF',int2str(DOF)))
%     hold on
%     demo = plot(0:1/120:1/120*(length(D)-1),D(:,DOF),'g-','LineWidth',4.5, 'DisplayName','demo');
%     smoothed = plot(0:1/120:1/120*(length(Dsmooth)-1),Dsmooth(:,DOF),'r-','LineWidth',3.5, 'DisplayName','35 smooth demo');
%     smoothed2 = plot(0:1/120:1/120*(length(Dsmooth2)-1),Dsmooth2(:,DOF),'b-','LineWidth',2.5, 'DisplayName','50 smooth demo demo');
%     smoothed3 = plot(0:1/120:1/120*(length(Dsmooth3)-1),Dsmooth3(:,DOF),'k-','LineWidth',1.5, 'DisplayName','65 smooth demo demo');
%     legend([demo, smoothed, smoothed2, smoothed3])
% end


%generate demo struct
demo_traj = generateDemo(q', 1/120);
%no sign flips as normal DMPs are exported before sign flips in DMP
%scripts!



%% Joint plots
% for joint = 1:7
%     figure('Name',strcat('joint vel',int2str(joint)))
%     hold on
%     plot_demo = plot(demo_traj.t,demo_traj.pos(joint,:),'g-','LineWidth',4.5, 'DisplayName','demo');
%     plot_TC_DMPpos = plot(0:1/1000:1/1000*(length(TC_DMPpos)-1),TC_DMPpos(:,joint),'b-','LineWidth',3.5, 'DisplayName','UC-TC-DMP');
%     plot_Opt_DMPpos = plot(0:1/1000:1/1000*(length(Opt_DMPpos)-1),Opt_DMPpos(:,joint),'r-','LineWidth',2.5, 'DisplayName','UC-Opt-DMP');
%     legend([plot_demo, plot_TC_DMPpos, plot_Opt_DMPpos])
% end


%% Joint vel plots
for joint = 1:7
    figure('Name',strcat('joint vel',int2str(joint)))
    hold on
    plot_demo = plot(demo_traj.t,demo_traj.vel(joint,:),'g-','LineWidth',4.5, 'DisplayName','demo');
    plot_TC_DMPvel = plot(0:1/1000:1/1000*(length(TC_DMPvel)-1),TC_DMPvel(:,joint),'b-','LineWidth',3.5, 'DisplayName','UC-TC-DMP');
    plot_Opt_DMPvel = plot(0:1/1000:1/1000*(length(Opt_DMPvel)-1),Opt_DMPvel(:,joint),'r-','LineWidth',2.5, 'DisplayName','UC-Opt-DMP');
    legend([plot_demo, plot_TC_DMPvel, plot_Opt_DMPvel])
end


%% Pose plots
% 
% poseDemo = ForwardKinematics(demo_traj.pos');
% poseTC_DMP = ForwardKinematics(TC_DMPpos);
% poseOpt_DMP = ForwardKinematics(Opt_DMPpos);
% 
% for DOF = 1:7
%     figure('Name',strcat('DOF',int2str(DOF)))
%     hold on
%     demo = plot(demo_traj.t,poseDemo(:,DOF),'g-','LineWidth',3.5, 'DisplayName','demo');
%     plot_TC_DMP = plot(0:1/1000:1/1000*(length(poseTC_DMP)-1),poseTC_DMP(:,DOF),'b-','LineWidth',2.5, 'DisplayName','TC_DMP');
%     plot_Opt_DMPvel = plot(0:1/1000:1/1000*(length(poseOpt_DMP)-1),poseOpt_DMP(:,DOF),'r-','LineWidth',1.5, 'DisplayName','Opt_DMP');
%     legend([demo, plot_TC_DMP, plot_Opt_DMPvel])
% end


