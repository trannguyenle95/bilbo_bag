clearvars
clc
close all 
addpath('../SupportScripts')
addpath('../Data/demos')
addpath('../Data/trajectories')

% Create demonstration trajectory
bag = 'A'
filename = '10l_bag_flip.csv';

UC_TC_DMPpos = csvread(strcat(bag, '_UC_TC_DMP_joint_10l_bag_flip.csv'));
UC_TC_DMPvel = csvread(strcat(bag, '_UC_TC_DMP_joint_vel_10l_bag_flip.csv'));
UC_Opt_DMPpos = csvread(strcat(bag, '_UC_Opt_joint_10l_bag_flip.csv'));
UC_Opt_DMPvel = csvread(strcat(bag, '_UC_Opt_joint_vel_10l_bag_flip.csv'));
%no sign flips as normal DMPs are exported before sign flips in DMP
%scripts!

generateJointDemo

%% Joint plots
figure('Name', 'Unconstrained DMP joint positions')
for joint = 1:7
    if joint < 5
        subplot(2,4,joint)
    else
        subplot(2,4,joint+0.5)
    end
    title(strcat('joint', num2str(joint), ' position'))
    hold on
    plot_demo = plot(demo_traj.t,demo_traj.pos(joint,:),'g-','LineWidth',7.5, 'DisplayName','demo');
    plot_UC_TC_DMPpos = plot(0:1/1000:1/1000*(length(UC_TC_DMPpos)-1),UC_TC_DMPpos(:,joint),'b-','LineWidth',4.5, 'DisplayName','UC-TC-DMP');
    plot_UC_Opt_DMPpos = plot(0:1/1000:1/1000*(length(UC_Opt_DMPpos)-1),UC_Opt_DMPpos(:,joint),'r-','LineWidth',1.5, 'DisplayName','UC-Opt-DMP');
    legend([plot_demo, plot_UC_TC_DMPpos, plot_UC_Opt_DMPpos], 'Location','northwest')
end


%% Joint vel plots
figure('Name', 'Unconstrained DMP joint velocities')
for joint = 1:7
    if joint < 5
        subplot(2,4,joint)
    else
        subplot(2,4,joint+0.5)
    end
    title(strcat('joint', num2str(joint), ' velocity'))
    hold on
    plot_demo = plot(demo_traj.t,demo_traj.vel(joint,:),'g-','LineWidth',7.5, 'DisplayName','demo');
    plot_UC_TC_DMPvel = plot(0:1/1000:1/1000*(length(UC_TC_DMPvel)-1),UC_TC_DMPvel(:,joint),'b-','LineWidth',4.5, 'DisplayName','UC-TC-DMP');
    plot_UC_Opt_DMPvel = plot(0:1/1000:1/1000*(length(UC_Opt_DMPvel)-1),UC_Opt_DMPvel(:,joint),'r-','LineWidth',1.5, 'DisplayName','UC-Opt-DMP');
    legend([plot_demo, plot_UC_TC_DMPvel, plot_UC_Opt_DMPvel], 'Location','northwest')
end


%% Pose plots
%>> TODO CHECK OUT WHETHER TO KEEP THIS CODE OR NOT
% poseDemo = ForwardKinematics(demo_traj.pos');
% poseTC_DMP = ForwardKinematics(UC_TC_DMPpos);
% poseOpt_DMP = ForwardKinematics(UC_Opt_DMPpos);
% 
% for DOF = 1:7
%     figure('Name',strcat('DOF',int2str(DOF)))
%     hold on
%     demo = plot(demo_traj.t,poseDemo(:,DOF),'g-','LineWidth',3.5, 'DisplayName','demo');
%     plot_TC_DMP = plot(0:1/1000:1/1000*(length(poseTC_DMP)-1),poseTC_DMP(:,DOF),'b-','LineWidth',2.5, 'DisplayName','TC_DMP');
%     plot_UC_Opt_DMPvel = plot(0:1/1000:1/1000*(length(poseOpt_DMP)-1),poseOpt_DMP(:,DOF),'r-','LineWidth',1.5, 'DisplayName','Opt_DMP');
%     legend([demo, plot_TC_DMP, plot_UC_Opt_DMPvel])
% end


%% Compare joint velocities and accelerations
v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
a_max = [10; 7.5; 10; 10; 10; 10; 10];

tau_DMPvel = csvread(strcat(bag, '_tau_DMP_joint_vel_10l_bag_flip.csv'))';
TC_DMPvel = csvread(strcat(bag, '_TC_DMP_joint_vel_10l_bag_flip.csv'))';
Opt_DMPvel = csvread(strcat(bag, '_Opt_DMP_joint_vel_10l_bag_flip.csv'))';

figure('Name', 'Joint velocities')
for joint = 1:7
    if joint < 5
        subplot(2,4,joint)
    else
        subplot(2,4,joint+0.5)
    end

    if mod(joint,2) == 1
        tau_DMPvel(joint,:) = -tau_DMPvel(joint,:);
        TC_DMPvel(joint,:) = -TC_DMPvel(joint,:);
        Opt_DMPvel(joint,:) = -Opt_DMPvel(joint,:);
    end

    title(strcat('joint', num2str(joint), ' velocity'))
    hold on
    plot_tau_DMP = plot(0:1/1000:1/1000*(length(tau_DMPvel)-1),tau_DMPvel(joint,:),'b-', 'LineWidth', 1.5, 'DisplayName','tau-DMP');
    plot_TC_DMP = plot(0:1/1000:1/1000*(length(TC_DMPvel)-1),TC_DMPvel(joint,:),'g-', 'LineWidth', 1.5, 'DisplayName','TC-DMP');
    plot_Opt_DMP = plot(0:1/1000:1/1000*(length(Opt_DMPvel)-1),Opt_DMPvel(joint,:),'k-', 'LineWidth', 1.5, 'DisplayName','Opt-DMP');

    limit_plot = plot([0 length(tau_DMPvel)/1000],[v_max(joint) v_max(joint)],'r:','LineWidth',1.5, 'DisplayName', 'Limit')
    plot([0 length(tau_DMPvel)/1000],-[v_max(joint) v_max(joint)],'r:','LineWidth',1.5)

    area([2,6],[v_max(joint),v_max(joint)], "FaceColor","y","FaceAlpha",0.3, "EdgeAlpha", 0)
    area([2,6],-[v_max(joint),v_max(joint)], "FaceColor","y","FaceAlpha",0.3, "EdgeAlpha", 0)

    legend([plot_tau_DMP, plot_TC_DMP, plot_Opt_DMP, limit_plot], 'Location','northwest')
end

tau_DMPacc = csvread(strcat(bag, '_tau_DMP_joint_acc_10l_bag_flip.csv'))';
TC_DMPacc = csvread(strcat(bag, '_TC_DMP_joint_acc_10l_bag_flip.csv'))';
Opt_DMPacc = csvread(strcat(bag, '_Opt_DMP_joint_acc_10l_bag_flip.csv'))';

figure('Name', 'Joint accelerations')
for joint = 1:7
    if joint < 5
        subplot(2,4,joint)
    else
        subplot(2,4,joint+0.5)
    end

    if mod(joint,2) == 1
        tau_DMPacc(joint,:) = -tau_DMPacc(joint,:);
        TC_DMPacc(joint,:) = -TC_DMPacc(joint,:);
        Opt_DMPacc(joint,:) = -Opt_DMPacc(joint,:);
    end

    title(strcat('joint', num2str(joint), ' acceleration'))
    hold on
    plot_tau_DMP = plot(0:1/1000:1/1000*(length(tau_DMPacc)-1),tau_DMPacc(joint,:),'b-', 'LineWidth', 1.5, 'DisplayName','tau-DMP');
    plot_TC_DMP = plot(0:1/1000:1/1000*(length(TC_DMPacc)-1),TC_DMPacc(joint,:),'g-', 'LineWidth', 1.5, 'DisplayName','TC-DMP');
    plot_Opt_DMP = plot(0:1/1000:1/1000*(length(Opt_DMPacc)-1),Opt_DMPacc(joint,:),'k-', 'LineWidth', 1.5, 'DisplayName','Opt-DMP');

    limit_plot = plot([0 length(tau_DMPacc)/1000],[v_max(joint) v_max(joint)],'r:','LineWidth',1.5, 'DisplayName', 'Limit')
    plot([0 length(tau_DMPacc)/1000],-[v_max(joint) v_max(joint)],'r:','LineWidth',1.5)

    area([2,6],[v_max(joint),v_max(joint)], "FaceColor","y","FaceAlpha",0.3, "EdgeAlpha", 0)
    area([2,6],-[v_max(joint),v_max(joint)], "FaceColor","y","FaceAlpha",0.3, "EdgeAlpha", 0)

    legend([plot_tau_DMP, plot_TC_DMP, plot_Opt_DMP, limit_plot], 'Location','northwest')
end
