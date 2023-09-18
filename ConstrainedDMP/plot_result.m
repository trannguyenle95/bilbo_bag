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

figure('Name','Manipulability measures')
subplot(5,1,1)
hold on
plot(0:1/120:1/120*(length(D)-1), D(:,1),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_x')
plot(0:1/120:1/120*(length(D)-1), D(:,2),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_y')
plot(0:1/120:1/120*(length(D)-1), D(:,3),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_z')
xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(5,1,2)
hold on
plot(0:1/120:1/120*(length(D)-1), mu1(:),'LineWidth',2.5, 'DisplayName','mu_1')
hold off
legend

subplot(5,1,3)
plot(0:1/120:1/120*(length(D)-1), mu2(:),'LineWidth',2.5, 'DisplayName','mu_2')
legend

subplot(5,1,4)
hold on
plot(0:1/120:1/120*(length(D)-1), mu3(:),'LineWidth',3.5, 'DisplayName','mu_3')
plot(0:1/120:1/120*(length(D)-1), w(:),'LineWidth',2, 'DisplayName','w')
hold off
legend

subplot(5,1,5)
hold on
plot(0:1/120:1/120*(length(D)-1), kappa(:),'LineWidth',2.5, 'DisplayName','kappa')
hold off
legend


%%%

figure('Name','Manipulability measures (traj after Constrained DMP)')
subplot(6,1,1)
hold on
plot(res{1}.t,poseDMP(:,1),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_x')
plot(res{1}.t,poseDMP(:,2),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_y')
plot(res{1}.t,poseDMP(:,3),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_z')
%xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(6,1,2)
hold on
plot(res{1}.t,res{1}.tau,'b','LineWidth',2.5, 'DisplayName', 'tau')
% plot(res{1}.t,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_{qx}')
% plot(res{1}.t,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_{qy}')
% plot(res{1}.t,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_{qz}')
% plot(res{1}.t,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMP_{qw}')
% xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off

subplot(6,1,3)
hold on
plot(res{1}.t, res_mu1(:),'LineWidth',2.5, 'DisplayName','mu_1, close to 1 better')
hold off
legend

subplot(6,1,4)
plot(res{1}.t, res_mu2(:),'LineWidth',2.5, 'DisplayName','mu_2, close to 1 better')
legend

subplot(6,1,5)
hold on
plot(res{1}.t, res_mu3(:),'LineWidth',3.5, 'DisplayName','mu_3, larger better')
plot(res{1}.t, res_w(:),'LineWidth',2, 'DisplayName','w')
hold off
legend

subplot(6,1,6)
hold on
plot(res{1}.t, res_kappa(:),'LineWidth',2.5, 'DisplayName','kappa, close to 1 better')
hold off
legend
