close all

%plots for all joints
for joint = 1:7
    plotJointTime(joint, res, sim_params, p_min, p_max, demo_traj)
end

% % Tau
% figure('Name','Tau')
% plot(res{1}.s,res{1}.tau,'b','LineWidth',2.5), hold on
% %plot(res{2}.s,res{2}.tau,'k--','LineWidth',2.5)
% %plot(res{3}.s,res{3}.tau,'m-.','LineWidth',2.5)
% %xticks([round(min([res{1}.s res{2}.s res{3}.s]),2) 0.9 1])
% xticks([round(min(res{1}.s),2) 0.9 1])
% yticks([])
% set(gca,'FontSize',16)
% %xlim([min([res{1}.s res{2}.s res{3}.s]) 1])
% xlim([min(res{1}.s) 1])
% ylim([0.9*dmp.nominal_tau 1.05*max(res{1}.tau)])
% xlabel('$s$', 'Interpreter','latex','Fontsize',30)
% ylabel('$\tau$', 'Interpreter','latex','Fontsize',30)
% set(gca,'XDir','reverse')

% Tau vs time (high tau means that slowed down motion is executed durnig these timesteps)
figure('Name','Tau')
plot(res{1}.t,res{1}.tau,'b','LineWidth',2.5), hold on
%xticks([round(min([res{1}.s res{2}.s res{3}.s]),2) 0.9 1])
%xticks([round(min(res{1}.s),2) 0.9 1])
%yticks([])
set(gca,'FontSize',16)
%xlim([min([res{1}.s res{2}.s res{3}.s]) 1])
%xlim([min(res{1}.s) 1])
ylim([0.9*dmp.nominal_tau 1.05*max(res{1}.tau)])
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
ylabel('$\tau$', 'Interpreter','latex','Fontsize',30)
%set(gca,'XDir','reverse')


%{
% s
figure('Name','s')
plot(res{1}.t,res{1}.s,'b','LineWidth',2.5), hold on
%plot(res{2}.t,res{2}.s,'k--','LineWidth',2.5)
%plot(res{3}.t,res{3}.s,'m-.','LineWidth',2.5)
%yticks([round(min([res{1}.s res{2}.s res{3}.s]),2) 0.9 1])
yticks([round(min(res{1}.s),2) 0.9 1])
set(gca,'FontSize',20)
%xlim([0 1.05*max([res{1}.t(end) res{2}.t(end) res{3}.t(end)])])
xlim([0 1.05*max(res{1}.t(end))])
ylim([0.9*res{1}.s(end) 1.05])
xlabel('$t$', 'Interpreter','latex','Fontsize',40)
ylabel('$s$', 'Interpreter','latex','Fontsize',40)
%}


figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

plot(res{1}.t,poseDMP(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_x')
plot(res{1}.t,poseDMP(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_y')
plot(res{1}.t,poseDMP(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(res{1}.t,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_{qx}')
plot(res{1}.t,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_{qy}')
plot(res{1}.t,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_{qz}')
plot(res{1}.t,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMP_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off


%%%


figure('Name','IKFK Cartesian Pose')
subplot(2,1,1)
hold on

plot(0:1/120:1/120*(length(D)-1),poseIKFK(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','IKFK_x')
plot(0:1/120:1/120*(length(D)-1),poseIKFK(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','IKFK_y')
plot(0:1/120:1/120*(length(D)-1),poseIKFK(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','IKFK_z')

plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1),poseIKFK(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','IKFK_{qx}')
plot(0:1/120:1/120*(length(D)-1),poseIKFK(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','IKFK_{qy}')
plot(0:1/120:1/120*(length(D)-1),poseIKFK(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','IKFK_{qz}')
plot(0:1/120:1/120*(length(D)-1),poseIKFK(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','IKFK_{qw}')

plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend

hold off

%%%


function plotJointTime(index, res, sim_params, p_min, p_max, demo_traj)
    figure('Name',strcat('joint',int2str(index)))
    sgtitle(strcat('joint',int2str(index)))
    subplot(3,1,1)
    hold on
    plot(res{1}.t,res{1}.ref_pos(index,:),'b-','LineWidth',1.5)
    plot(demo_traj.t,demo_traj.pos(index,:),'r-','LineWidth',1.5)

    plot([res{1}.t(1) res{1}.t(end)],[p_max(index) p_max(index)],'r:','LineWidth',1)
    plot([res{1}.t(1) res{1}.t(end)],[p_min(index) p_min(index)],'r:','LineWidth',1)
    %yticks(sort([p_min(index) 0 p_max(index)]))
    %TODO: auto ticks or specified? both?
    
    %xticks([round(min(res{1}.t),2) 0.9 1])
    %xlim([min(res{1}.t) 1])
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    %ylim([min(res{1}.ref_pos(index,:)) max(res{1}.ref_pos(index,:))])
    ylim([min([res{1}.ref_pos(index,:) demo_traj.pos(index,:)]) max([res{1}.ref_pos(index,:) demo_traj.pos(index,:)])])

    ylabel('$y$', 'Interpreter','latex','Fontsize',30)
    %set(gca,'XDir','reverse')
    set(gca,'FontSize',16)
    hold off
    
    subplot(3,1,2)
    hold on
    plot(res{1}.t,res{1}.ref_vel(index,:),'b-','LineWidth',1.5)
    plot(demo_traj.t,demo_traj.vel(index,:),'r-','LineWidth',1.5)
    plot([res{1}.t(1) res{1}.t(end)],[sim_params.v_max(index) sim_params.v_max(index)],'r:','LineWidth',1)
    plot([res{1}.t(1) res{1}.t(end)],-[sim_params.v_max(index) sim_params.v_max(index)],'r:','LineWidth',1)
    %yticks([-sim_params.v_max(index) 0 sim_params.v_max(index)])
    
    %xticks([round(min(res{1}.t),2) 0.9 1])
    %xlim([min(res{1}.t) 1])
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min(res{1}.ref_vel(index,:)) max(res{1}.ref_vel(index,:))])
    %ylim([min([res{1}.ref_vel(index,:) demo_traj.vel(index,:)]) max([res{1}.ref_vel(index,:) demo_traj.vel(index,:)])])

    ylabel('$\dot{y}$', 'Interpreter','latex','Fontsize',30)
    %set(gca,'XDir','reverse')
    set(gca,'FontSize',16)
    hold off
    
    subplot(3,1,3)
    hold on
    plot(res{1}.t,res{1}.ref_acc(index,:),'b-','LineWidth',1.5)
    plot(demo_traj.t,demo_traj.acc(index,:),'r-','LineWidth',1.5)
    plot([res{1}.t(1) res{1}.t(end)],[sim_params.a_max(index) sim_params.a_max(index)],'r:','LineWidth',1)
    plot([res{1}.t(1) res{1}.t(end)],-[sim_params.a_max(index) sim_params.a_max(index)],'r:','LineWidth',1)
    %yticks([-sim_params.a_max(index) 0 sim_params.a_max(index)])
    
    %xticks([round(min(res{1}.t),2) 0.9 1])
    %xlim([min(res{1}.t) 1])
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min(res{1}.ref_acc(index,:)) max(res{1}.ref_acc(index,:))])
    %ylim([min([res{1}.ref_acc(index,:) demo_traj.acc(index,:)]) max([res{1}.ref_acc(index,:) demo_traj.acc(index,:)])])

    ylabel('$\ddot{y}$', 'Interpreter','latex','Fontsize',30)
    %set(gca,'XDir','reverse')
    set(gca,'FontSize',16)
    hold off
end