%% Joint plots

%plots for all joints
for joint = 1:7
    plotJointTime(joint, data{version}, actual_pos_lim , actual_vel_lim, actual_accel_lim, data{2})
end


%% Cartesian plots
%index 4 = pos optimization
figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(data{2}.Time,poseUnconstrainedDMP(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','UC-DMP_x')
plot(data{2}.Time,poseUnconstrainedDMP(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','UC-DMP_y')
plot(data{2}.Time,poseUnconstrainedDMP(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','UC-DMP_z')


plot(data{4}.Time,poseDMP(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','Opt-DMP_x')
plot(data{4}.Time,poseDMP(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','Opt-DMP_y')
plot(data{4}.Time,poseDMP(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','Opt-DMP_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(data{2}.Time,poseUnconstrainedDMP(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','UC-DMP_{qx}')
plot(data{2}.Time,poseUnconstrainedDMP(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','UC-DMP_{qy}')
plot(data{2}.Time,poseUnconstrainedDMP(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','UC-DMP_{qz}')
plot(data{2}.Time,poseUnconstrainedDMP(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','UC-DMP_{qw}')

plot(data{4}.Time,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','Opt-DMP_{qx}')
plot(data{4}.Time,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','Opt-DMP_{qy}')
plot(data{4}.Time,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','Opt-DMP_{qz}')
plot(data{4}.Time,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','Opt-DMP_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

%% Plot DMP, both constrained DMPs and Demo in same plot, like in the paper code comes from
% label_font = 17;
% ax_fontsize = 14;
% ind = [1 2 3 4 5 6 7]; % choose DoFs to plot
% for k=1:length(ind)
%     i = ind(k);
%     fig = figure;
%     sgtitle(strcat('joint',int2str(i)))
% 
%     % plot joint positions
%     ax = subplot(3,1,1);
%     hold on;
%     legend_ = {};
%     for k=1:length(data)
%         if (~data{k}.plot2D), continue; end
%         plot(data{k}.Time, data{k}.Pos(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
%         legend_ = [legend_ data{k}.legend];
%     end
%     axis tight;
%     plot(ax.XLim, [actual_pos_lim(i,1) actual_pos_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
%     plot(ax.XLim, [actual_pos_lim(i,2) actual_pos_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
%     ylabel('pos [$rad$]', 'interpreter','latex', 'fontsize',label_font);
%     legend(legend_, 'interpreter','latex', 'fontsize',17, 'Position',[0.2330 0.9345 0.5520 0.0294], 'Orientation', 'horizontal');
%     ax.FontSize = ax_fontsize;
%     hold off;
% 
%     % plot joint velocities
%     ax = subplot(3,1,2);
%     hold on;
%     for k=1:length(data)
%         if (~data{k}.plot2D), continue; end
%         plot(data{k}.Time, data{k}.Vel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
%     end
%     axis tight;
%     plot(ax.XLim, [actual_vel_lim(i,1) actual_vel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
%     plot(ax.XLim, [actual_vel_lim(i,2) actual_vel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
%     ylabel('vel [$rad/s$]', 'interpreter','latex', 'fontsize',label_font);
%     ax.FontSize = ax_fontsize;
%     hold off;
% 
%     % plot joint velocities
%     ax = subplot(3,1,3);
%     hold on;
%     for k=1:length(data)
%         if (~data{k}.plot2D), continue; end
%         plot(data{k}.Time, data{k}.Accel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
%     end
%     axis tight;
%     plot(ax.XLim, [actual_accel_lim(i,1) actual_accel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
%     plot(ax.XLim, [actual_accel_lim(i,2) actual_accel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
%     ylabel('accel [$rad/s^2$]', 'interpreter','latex', 'fontsize',label_font);
%     xlabel('time [$s$]', 'interpreter','latex', 'fontsize',label_font);
%     ax.FontSize = ax_fontsize;
%     hold off;
% end

% ======================================================


%% 
function plotJointTime(joint, traj, p_lim , v_lim, a_lim, unconstrained_DMP)
    figure('Name',strcat('joint',int2str(joint)))
    sgtitle(strcat('joint',int2str(joint)))

    subplot(3,1,1)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.Time, unconstrained_DMP.Pos(joint,:), 'r-', 'LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(traj.Time, traj.Pos(joint,:), 'b-', 'LineWidth',1.5, 'DisplayName','Opt-DMP');
    plot([traj.Time(1) traj.Time(end)],[p_lim(joint,1) p_lim(joint,1)],'r:','LineWidth',1) %lower limit
    plot([traj.Time(1) traj.Time(end)],[p_lim(joint,2) p_lim(joint,2)],'r:','LineWidth',1) %upper limit
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([traj.Pos(joint,:) unconstrained_DMP.Pos(joint,:) p_lim(joint,1)]) max([traj.Pos(joint,:) unconstrained_DMP.Pos(joint,:) p_lim(joint,2)])])
    xlim([0 traj.Time(end)])
    ylabel('$y$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off


    subplot(3,1,2)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.Time, unconstrained_DMP.Vel(joint,:), 'r-', 'LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(traj.Time, traj.Vel(joint,:), 'b-', 'LineWidth',1.5, 'DisplayName','Opt-DMP');   
    plot([traj.Time(1) traj.Time(end)],[v_lim(joint,1) v_lim(joint,1)],'r:','LineWidth',1) %lower limit
    plot([traj.Time(1) traj.Time(end)],[v_lim(joint,2) v_lim(joint,2)],'r:','LineWidth',1) %upper limit
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([traj.Vel(joint,:) unconstrained_DMP.Vel(joint,:) v_lim(joint,1)]) max([traj.Vel(joint,:) unconstrained_DMP.Vel(joint,:) v_lim(joint,2)])])
    xlim([0 traj.Time(end)])
    ylabel('$\dot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off


    subplot(3,1,3)
    hold on
    unconstrained_dmp_plot = plot(unconstrained_DMP.Time, unconstrained_DMP.Accel(joint,:), 'r-', 'LineWidth',1.5, 'DisplayName','UC-DMP');
    constrained_dmp_plot = plot(traj.Time, traj.Accel(joint,:), 'b-', 'LineWidth',1.5, 'DisplayName','Opt-DMP');   
    plot([traj.Time(1) traj.Time(end)],[a_lim(joint,1) a_lim(joint,1)],'r:','LineWidth',1) %lower limit
    plot([traj.Time(1) traj.Time(end)],[a_lim(joint,2) a_lim(joint,2)],'r:','LineWidth',1) %upper limit
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([traj.Accel(joint,:) unconstrained_DMP.Accel(joint,:) a_lim(joint,1)]) max([traj.Accel(joint,:) unconstrained_DMP.Accel(joint,:) a_lim(joint,2)])])
    xlim([0 traj.Time(end)])
    ylabel('$\ddot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([unconstrained_dmp_plot, constrained_dmp_plot])
    hold off

end