close all

%DEFINE FUNCTION FOR PLOTTING (joint index as input param!)
for joint = 1:7
    plotJointTime(joint, res, sim_params, p_min, p_max, demo_traj)
end

% Tau
figure('Name','Tau')
plot(res{1}.s,res{1}.tau,'b','LineWidth',2.5), hold on
%plot(res{2}.s,res{2}.tau,'k--','LineWidth',2.5)
%plot(res{3}.s,res{3}.tau,'m-.','LineWidth',2.5)
%xticks([round(min([res{1}.s res{2}.s res{3}.s]),2) 0.9 1])
xticks([round(min(res{1}.s),2) 0.9 1])
yticks([])
set(gca,'FontSize',16)
%xlim([min([res{1}.s res{2}.s res{3}.s]) 1])
xlim([min(res{1}.s) 1])
ylim([0.9*dmp.nominal_tau 1.05*max(res{1}.tau)])
xlabel('$s$', 'Interpreter','latex','Fontsize',30)
ylabel('$\tau$', 'Interpreter','latex','Fontsize',30)
set(gca,'XDir','reverse')

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
%multiply quaternion by -1 to get same representation as in DMP for
%comparison
plot(0:1/120:1/120*(length(D)-1), -D(:,4),'color','#DFE916','LineWidth',2.5, 'DisplayName','demo_{qx}')
plot(0:1/120:1/120*(length(D)-1), -D(:,5),'color','#00FF1B','LineWidth',2.5, 'DisplayName','demo_{qy}')
plot(0:1/120:1/120*(length(D)-1), -D(:,6),'color','#49A9B6','LineWidth',2.5, 'DisplayName','demo_{qz}')
plot(0:1/120:1/120*(length(D)-1), -D(:,7),'color','#7A364F','LineWidth',2.5, 'DisplayName','demo_{qw}')

plot(res{1}.t,poseDMP(:,4),'color','#2016E9','LineWidth',1.5, 'DisplayName','DMP_{qx}')
plot(res{1}.t,poseDMP(:,5),'color','#FF00E4','LineWidth',1.5, 'DisplayName','DMP_{qy}')
plot(res{1}.t,poseDMP(:,6),'color','#B65649','LineWidth',1.5, 'DisplayName','DMP_{qz}')
plot(res{1}.t,poseDMP(:,7),'color','#367A62','LineWidth',1.5, 'DisplayName','DMP_{qw}')

xlabel('$t$', 'Interpreter','latex','Fontsize',30)
legend
hold off
