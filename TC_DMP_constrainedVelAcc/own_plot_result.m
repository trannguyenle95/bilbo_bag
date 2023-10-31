close all
%plots for all joints
for joint = 1:7
    plotJointTime(joint, res, sim_params, p_min, p_max, demo_traj)
end

%%%
% Tau vs time (high tau means that slowed down motion is executed durnig these timesteps)
figure('Name','Tau')
plot(res.t,res.tau,'b','LineWidth',2.5), hold on
set(gca,'FontSize',16)
ylim([0.9*dmp.nominal_tau 1.05*max(res.tau)])
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
ylabel('$\tau$', 'Interpreter','latex','Fontsize',30)

%%%
figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

plot(res.t,poseDMP(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_x')
plot(res.t,poseDMP(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_y')
plot(res.t,poseDMP(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(res.t,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_{qx}')
plot(res.t,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_{qy}')
plot(res.t,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_{qz}')
plot(res.t,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMP_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off


%%%
function plotJointTime(index, res, sim_params, p_min, p_max, demo_traj)
    figure('Name',strcat('joint',int2str(index)))
    sgtitle(strcat('joint',int2str(index)))

    subplot(3,1,1)
    hold on
    dmp = plot(res.t,res.ref_pos(index,:),'b-','LineWidth',1.5, 'DisplayName','TC-DMP');
    demo = plot(demo_traj.t,demo_traj.pos(index,:),'r-','LineWidth',1.5, 'DisplayName','demo');
    plot([res.t(1) res.t(end)],[p_max(index) p_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],[p_min(index) p_min(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.ref_pos(index,:) demo_traj.pos(index,:) p_min(index)]) max([res.ref_pos(index,:) demo_traj.pos(index,:) p_max(index)])])
    xlim([0 max([res.t demo_traj.t])])
    ylabel('$y$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([dmp, demo])
    hold off
    
    subplot(3,1,2)
    hold on
    dmp = plot(res.t,res.ref_vel(index,:),'b-','LineWidth',1.5, 'DisplayName','TC-DMP');
    demo = plot(demo_traj.t,demo_traj.vel(index,:),'r-','LineWidth',1.5, 'DisplayName','demo');
    plot([res.t(1) res.t(end)],[sim_params.v_max(index) sim_params.v_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],-[sim_params.v_max(index) sim_params.v_max(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.ref_vel(index,:) demo_traj.vel(index,:) -sim_params.v_max(index)]) max([res.ref_vel(index,:) demo_traj.vel(index,:) sim_params.v_max(index)])])
    xlim([0 max([res.t demo_traj.t])])
    ylabel('$\dot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([dmp, demo])
    hold off
    
    subplot(3,1,3)
    hold on
    dmp = plot(res.t,res.ref_acc(index,:),'b-','LineWidth',1.5, 'DisplayName','TC-DMP');
    demo = plot(demo_traj.t,demo_traj.acc(index,:),'r-','LineWidth',1.5, 'DisplayName','demo');
    plot([res.t(1) res.t(end)],[sim_params.a_max(index) sim_params.a_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],-[sim_params.a_max(index) sim_params.a_max(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.ref_acc(index,:) demo_traj.acc(index,:) -sim_params.a_max(index)]) max([res.ref_acc(index,:) demo_traj.acc(index,:) sim_params.a_max(index)])])
    xlim([0 max([res.t demo_traj.t])])
    ylabel('$\ddot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([dmp, demo])
    hold off
end