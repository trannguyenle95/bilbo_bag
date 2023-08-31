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