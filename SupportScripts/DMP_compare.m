clearvars
clc
close all 

addpath('/home/erichannus/catkin_ws/src/Data/trajectories/')
DMP2pos = csvread('nominalDMP2_joint_10l_bag_flip.csv');
DMP2vel = csvread('nominalDMP2_joint_vel_10l_bag_flip.csv');
DMP3pos = csvread('nominalDMP3_joint_10l_bag_flip.csv');
DMP3vel = csvread('nominalDMP3_joint_vel_10l_bag_flip.csv');


addpath('/home/erichannus/catkin_ws/src/SupportScripts/')
addpath('/home/erichannus/catkin_ws/src/Data/demos/')

% Create demonstration trajectory
filename = '10l_bag_flip.csv';
D = preprocess(filename, false, 0.00, 0.00, 0.00, 1, 'ori1', 0.38);
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


%test that IK produces same result
% for i = 1:10
%     [q, jacobians] = InverseKinematics(Dsmooth);
%     demo_traj = generateDemo(q', 1/120);
%     figure('Name',strcat('iter',int2str(i)))
%     plot(demo_traj.t,demo_traj.pos(2,:),'g-','LineWidth',3.5, 'DisplayName','demo');
% end


%% Joint plots
% for joint = 1:7
%     figure('Name',strcat('joint',int2str(joint)))
%     hold on
%     demo = plot(demo_traj.t,demo_traj.pos(joint,:),'g-','LineWidth',4.5, 'DisplayName','demo');
%     dmp2 = plot(0:1/1000:1/1000*(length(DMP2pos)-1),DMP2pos(:,joint),'b-','LineWidth',3.5, 'DisplayName','DMP2');
%     dmp3 = plot(0:1/1000:1/1000*(length(DMP3pos)-1),DMP3pos(:,joint),'r-','LineWidth',2.5, 'DisplayName','DMP3');
%     legend([demo, dmp2, dmp3])
% end


%% Joint vel plots
for joint = 1:7
    figure('Name',strcat('joint vel',int2str(joint)))
    hold on
    demo = plot(demo_traj.t,demo_traj.vel(joint,:),'g-','LineWidth',4.5, 'DisplayName','demo');
    dmp2 = plot(0:1/1000:1/1000*(length(DMP2vel)-1),DMP2vel(:,joint),'b-','LineWidth',3.5, 'DisplayName','DMP2');
    dmp3 = plot(0:1/1000:1/1000*(length(DMP3vel)-1),DMP3vel(:,joint),'r-','LineWidth',2.5, 'DisplayName','DMP3');
    legend([demo, dmp2, dmp3])
end

%% Pose plots
% 
% poseDemo = ForwardKinematics(demo_traj.pos');
% poseDMP2 = ForwardKinematics(DMP2pos);
% poseDMP3 = ForwardKinematics(DMP3pos);
% 
% 
% 
% for DOF = 1:7
%     figure('Name',strcat('DOF',int2str(DOF)))
%     hold on
%     demo = plot(demo_traj.t,poseDemo(:,DOF),'g-','LineWidth',3.5, 'DisplayName','demo');
%     dmp2 = plot(0:1/1000:1/1000*(length(poseDMP2)-1),poseDMP2(:,DOF),'b-','LineWidth',2.5, 'DisplayName','DMP2');
%     dmp3 = plot(0:1/1000:1/1000*(length(poseDMP3)-1),poseDMP3(:,DOF),'r-','LineWidth',1.5, 'DisplayName','DMP3');
%     legend([demo, dmp2, dmp3])
% end


