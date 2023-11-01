close all
%plots for all joints
for joint = 1:7
    plotJointTime(joint, res, sim_params, p_min, p_max, nominalTraj)
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
plot(nominalTraj.t,poseUnconstrainedDMP(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','UC-DMP_x')
plot(nominalTraj.t,poseUnconstrainedDMP(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','UC-DMP_y')
plot(nominalTraj.t,poseUnconstrainedDMP(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','UC-DMP_z')
plot(res.t,poseDMP(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','TC-DMP_x')
plot(res.t,poseDMP(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','TC-DMP_y')
plot(res.t,poseDMP(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','TC-DMP_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(nominalTraj.t,poseUnconstrainedDMP(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','UC-DMP_{qx}')
plot(nominalTraj.t,poseUnconstrainedDMP(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','UC-DMP_{qy}')
plot(nominalTraj.t,poseUnconstrainedDMP(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','UC-DMP_{qz}')
plot(nominalTraj.t,poseUnconstrainedDMP(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','UC-DMP_{qw}')
plot(res.t,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','TC-DMP_{qx}')
plot(res.t,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','TC-DMP_{qy}')
plot(res.t,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','TC-DMP_{qz}')
plot(res.t,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','TC-DMP_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off


%%%
function plotJointTime(index, res, sim_params, p_min, p_max, unconstrained_DMP)
    figure('Name',strcat('joint',int2str(index)))
    sgtitle(strcat('joint',int2str(index)))

    subplot(3,1,1)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.t,unconstrained_DMP.pos(index,:),'r-','LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(res.t,res.ref_pos(index,:),'b-','LineWidth',1.5, 'DisplayName','TC-DMP');
    plot([res.t(1) res.t(end)],[p_max(index) p_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],[p_min(index) p_min(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.ref_pos(index,:) unconstrained_DMP.pos(index,:) p_min(index)]) max([res.ref_pos(index,:) unconstrained_DMP.pos(index,:) p_max(index)])])
    xlim([0 max([res.t unconstrained_DMP.t])])
    ylabel('$y$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off
    
    subplot(3,1,2)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.t,unconstrained_DMP.vel(index,:),'r-','LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(res.t,res.ref_vel(index,:),'b-','LineWidth',1.5, 'DisplayName','TC-DMP');
    plot([res.t(1) res.t(end)],[sim_params.v_max(index) sim_params.v_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],-[sim_params.v_max(index) sim_params.v_max(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.ref_vel(index,:) unconstrained_DMP.vel(index,:) -sim_params.v_max(index)]) max([res.ref_vel(index,:) unconstrained_DMP.vel(index,:) sim_params.v_max(index)])])
    xlim([0 max([res.t unconstrained_DMP.t])])
    ylabel('$\dot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off
    
    subplot(3,1,3)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.t,unconstrained_DMP.acc(index,:),'r-','LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(res.t,res.ref_acc(index,:),'b-','LineWidth',1.5, 'DisplayName','TC-DMP');
    plot([res.t(1) res.t(end)],[sim_params.a_max(index) sim_params.a_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],-[sim_params.a_max(index) sim_params.a_max(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.ref_acc(index,:) unconstrained_DMP.acc(index,:) -sim_params.a_max(index)]) max([res.ref_acc(index,:) unconstrained_DMP.acc(index,:) sim_params.a_max(index)])])
    xlim([0 max([res.t unconstrained_DMP.t])])
    ylabel('$\ddot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off
end