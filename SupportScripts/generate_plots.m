close all
%plots for all joints
for joint = 1:7
    plotJointTime(joint, res, p_min, p_max, v_max, a_max, unconstrained_DMP)
end


%%%
figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(unconstrained_DMP.t,poseUnconstrainedDMP(:,1), 'r-','LineWidth',2.5, 'DisplayName','UC-DMP_x')
plot(unconstrained_DMP.t,poseUnconstrainedDMP(:,2),'g-','LineWidth',2.5, 'DisplayName','UC-DMP_y')
plot(unconstrained_DMP.t,poseUnconstrainedDMP(:,3),'b-','LineWidth',2.5, 'DisplayName','UC-DMP_z')
plot(res.t,poseDMP(:,1),'r:','LineWidth',2.5, 'DisplayName', strcat(res.version, '_x'))
plot(res.t,poseDMP(:,2),'g:','LineWidth',2.5, 'DisplayName', strcat(res.version, '_y'))
plot(res.t,poseDMP(:,3),'b:','LineWidth',2.5, 'DisplayName', strcat(res.version, '_z'))
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(unconstrained_DMP.t,poseUnconstrainedDMP(:,4),'r-','LineWidth',2.5, 'DisplayName','UC-DMP_{qx}')
plot(unconstrained_DMP.t,poseUnconstrainedDMP(:,5),'g-','LineWidth',2.5, 'DisplayName','UC-DMP_{qy}')
plot(unconstrained_DMP.t,poseUnconstrainedDMP(:,6),'b-','LineWidth',2.5, 'DisplayName','UC-DMP_{qz}')
plot(unconstrained_DMP.t,poseUnconstrainedDMP(:,7),'k-','LineWidth',2.5, 'DisplayName','UC-DMP_{qv}')
plot(res.t,poseDMP(:,4),'r:','LineWidth',2.5, 'DisplayName', strcat(res.version, '_{qx}'))
plot(res.t,poseDMP(:,5),'g:','LineWidth',2.5, 'DisplayName',strcat(res.version, '_{qy}'))
plot(res.t,poseDMP(:,6),'b:','LineWidth',2.5, 'DisplayName',strcat(res.version, '_{qz}'))
plot(res.t,poseDMP(:,7),'k:','LineWidth',2.5, 'DisplayName',strcat(res.version, '_{qv}'))

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

%% Joint Velocity plots
figure('Name', 'Joint velocities')

for joint = 1:7
    if joint < 5
        subplot(2,4,joint)
    else
        subplot(2,4,joint+0.5)
    end
    title(strcat('joint', num2str(joint), ' velocity'))
    hold on

    plot_unconstrained = plot(unconstrained_DMP.t,unconstrained_DMP.vel(joint,:),'r-','LineWidth',1.5, 'DisplayName','UC-DMP');
    plot_tau_DMP = plot(res.t,res.vel(joint,:),'b-','LineWidth',1.5, 'DisplayName',res.version);
    limit_plot = plot([res.t(1) res.t(end)],[v_max(joint) v_max(joint)],'r:','LineWidth',1.5, 'DisplayName', 'Limit')
    plot([res.t(1) res.t(end)],-[v_max(joint) v_max(joint)],'r:','LineWidth',1.5)

    legend([plot_unconstrained, plot_tau_DMP, limit_plot], 'Location','northwest')
    %legend([plot_tau_DMP, limit_plot], 'Location','northwest')
end


%% Joint Acceleration plots
figure('Name', 'Joint accelerations')

for joint = 1:7
    if joint < 5
        subplot(2,4,joint)
    else
        subplot(2,4,joint+0.5)
    end
    title(strcat('joint', num2str(joint), ' acceleration'))
    hold on

    plot_unconstrained = plot(unconstrained_DMP.t,unconstrained_DMP.acc(joint,:),'r-','LineWidth',1.5, 'DisplayName','UC-DMP');
    plot_tau_DMP = plot(res.t,res.acc(joint,:),'b-','LineWidth',1.5, 'DisplayName',res.version);
    limit_plot = plot([res.t(1) res.t(end)],[a_max(joint) a_max(joint)],'r:','LineWidth',1.5, 'DisplayName', 'Limit')
    plot([res.t(1) res.t(end)],-[a_max(joint) a_max(joint)],'r:','LineWidth',1.5)

    legend([plot_unconstrained, plot_tau_DMP, limit_plot], 'Location','northwest')
    %legend([plot_tau_DMP, limit_plot], 'Location','northwest')
end


%%
function plotJointTime(index, res, p_min, p_max, v_max, a_max, unconstrained_DMP)
    figure('Name',strcat('joint',int2str(index)))
    sgtitle(strcat('joint',int2str(index)))

    subplot(3,1,1)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.t,unconstrained_DMP.pos(index,:),'r-','LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(res.t,res.pos(index,:),'b-','LineWidth',1.5, 'DisplayName',res.version);
    plot([res.t(1) res.t(end)],[p_max(index) p_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],[p_min(index) p_min(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.pos(index,:) unconstrained_DMP.pos(index,:) p_min(index)]) max([res.pos(index,:) unconstrained_DMP.pos(index,:) p_max(index)])])
    xlim([0 max([res.t unconstrained_DMP.t])])
    ylabel('$y$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off
    
    subplot(3,1,2)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.t,unconstrained_DMP.vel(index,:),'r-','LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(res.t,res.vel(index,:),'b-','LineWidth',1.5, 'DisplayName',res.version);
    plot([res.t(1) res.t(end)],[v_max(index) v_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],-[v_max(index) v_max(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.vel(index,:) unconstrained_DMP.vel(index,:) -v_max(index)]) max([res.vel(index,:) unconstrained_DMP.vel(index,:) v_max(index)])])
    xlim([0 max([res.t unconstrained_DMP.t])])
    ylabel('$\dot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off
    
    subplot(3,1,3)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.t,unconstrained_DMP.acc(index,:),'r-','LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(res.t,res.acc(index,:),'b-','LineWidth',1.5, 'DisplayName',res.version);
    plot([res.t(1) res.t(end)],[a_max(index) a_max(index)],'r:','LineWidth',1)
    plot([res.t(1) res.t(end)],-[a_max(index) a_max(index)],'r:','LineWidth',1)
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([res.acc(index,:) unconstrained_DMP.acc(index,:) -a_max(index)]) max([res.acc(index,:) unconstrained_DMP.acc(index,:) a_max(index)])])
    xlim([0 max([res.t unconstrained_DMP.t])])
    ylabel('$\ddot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off
end