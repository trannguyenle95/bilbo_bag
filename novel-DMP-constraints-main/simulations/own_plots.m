%% Joint plots

%plots for all joints
for joint = 1:7
    plotJointTime(joint, data{version}, actual_pos_lim , actual_vel_lim, actual_accel_lim, demo_traj)
end


%% Cartesian plots


%index 2 = DMP without constraints
figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

plot(data{2}.Time,poseDMP2(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_x')
plot(data{2}.Time,poseDMP2(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_y')
plot(data{2}.Time,poseDMP2(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(data{2}.Time,poseDMP2(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_{qx}')
plot(data{2}.Time,poseDMP2(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_{qy}')
plot(data{2}.Time,poseDMP2(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_{qz}')
plot(data{2}.Time,poseDMP2(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMP_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

%index 3 = vel optimization
figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

plot(data{3}.Time,poseDMP3(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMPv_x')
plot(data{3}.Time,poseDMP3(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMPv_y')
plot(data{3}.Time,poseDMP3(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMPv_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(data{3}.Time,poseDMP3(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMPv_{qx}')
plot(data{3}.Time,poseDMP3(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMPv_{qy}')
plot(data{3}.Time,poseDMP3(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMPv_{qz}')
plot(data{3}.Time,poseDMP3(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMPv_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

%index 4 = pos optimization

figure('Name','Cartesian Pose')
subplot(2,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')

plot(data{4}.Time,poseDMP4(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMPp_x')
plot(data{4}.Time,poseDMP4(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMPp_y')
plot(data{4}.Time,poseDMP4(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMPp_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(2,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(data{4}.Time,poseDMP4(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMPp_{qx}')
plot(data{4}.Time,poseDMP4(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMPp_{qy}')
plot(data{4}.Time,poseDMP4(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMPp_{qz}')
plot(data{4}.Time,poseDMP4(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMPp_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

%% 
function plotJointTime(joint, traj, p_lim , v_lim, a_lim, demo)
    figure('Name',strcat('joint',int2str(joint)))
    sgtitle(strcat('joint',int2str(joint)))

    subplot(3,1,1)
    hold on

    dmp_plot = plot(traj.Time, traj.Pos(joint,:), 'b-', 'LineWidth',1.5, 'DisplayName','OptDMP');
    demo_plot = plot(demo.t,demo.pos(joint,:),'r-','LineWidth',1.5, 'DisplayName','demo');
       

    plot([traj.Time(1) traj.Time(end)],[p_lim(joint,1) p_lim(joint,1)],'r:','LineWidth',1) %lower limit
    plot([traj.Time(1) traj.Time(end)],[p_lim(joint,2) p_lim(joint,2)],'r:','LineWidth',1) %upper limit
    
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([traj.Pos(joint,:) demo.pos(joint,:) p_lim(joint,1)]) max([traj.Pos(joint,:) demo.pos(joint,:) p_lim(joint,2)])])
    xlim([0 max([traj.Time demo.t])])
    ylabel('$y$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([dmp_plot, demo_plot])
    hold off


    subplot(3,1,2)
    hold on

    dmp_plot = plot(traj.Time, traj.Vel(joint,:), 'b-', 'LineWidth',1.5, 'DisplayName','OptDMP');
    demo_plot = plot(demo.t,demo.vel(joint,:),'r-','LineWidth',1.5, 'DisplayName','demo');
       

    plot([traj.Time(1) traj.Time(end)],[v_lim(joint,1) v_lim(joint,1)],'r:','LineWidth',1) %lower limit
    plot([traj.Time(1) traj.Time(end)],[v_lim(joint,2) v_lim(joint,2)],'r:','LineWidth',1) %upper limit
    
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([traj.Vel(joint,:) demo.vel(joint,:) v_lim(joint,1)]) max([traj.Vel(joint,:) demo.vel(joint,:) v_lim(joint,2)])])
    xlim([0 max([traj.Time demo.t])])
    ylabel('$\dot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([dmp_plot, demo_plot])
    hold off
    

    subplot(3,1,3)
    hold on

    dmp_plot = plot(traj.Time, traj.Accel(joint,:), 'b-', 'LineWidth',1.5, 'DisplayName','OptDMP');
    demo_plot = plot(demo.t,demo.acc(joint,:),'r-','LineWidth',1.5, 'DisplayName','demo');
       

    plot([traj.Time(1) traj.Time(end)],[a_lim(joint,1) a_lim(joint,1)],'r:','LineWidth',1) %lower limit
    plot([traj.Time(1) traj.Time(end)],[a_lim(joint,2) a_lim(joint,2)],'r:','LineWidth',1) %upper limit
    
    xlabel('$t$', 'Interpreter','latex','Fontsize',30)
    ylim([min([traj.Accel(joint,:) demo.acc(joint,:) a_lim(joint,1)]) max([traj.Accel(joint,:) demo.acc(joint,:) a_lim(joint,2)])])
    xlim([0 max([traj.Time demo.t])])
    ylabel('$\ddot{y}$', 'Interpreter','latex','Fontsize',30)
    set(gca,'FontSize',16)
    legend([dmp_plot, demo_plot])
    hold off

end