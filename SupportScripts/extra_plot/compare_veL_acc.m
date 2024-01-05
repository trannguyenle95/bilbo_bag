clear all
clc

%the .mat files were all saved from the trajectories fitted to bag A
unconstrained_traj = load("unconstrained_traj.mat")
Opt_DMP_traj = load("Opt_DMP_traj.mat")
%Opt_DMP_traj = load("vel_constrained_Opt_DMP.mat") %ONLY USE THIS FOR COMPARISON - not actually used in experiments!

tau_DMP_res = load("tau_DMP_res.mat")
TC_DMP_res = load("TC_DMP_res.mat")

%make uniform structs for easier plotting
unconstrained.t = unconstrained_traj.unconstrained.Time
unconstrained.p = unconstrained_traj.unconstrained.Pos
unconstrained.v = unconstrained_traj.unconstrained.Vel
unconstrained.a = unconstrained_traj.unconstrained.Accel

tau_DMP.t = tau_DMP_res.res.t
tau_DMP.p = tau_DMP_res.res.pos
tau_DMP.v = tau_DMP_res.res.vel
tau_DMP.a = tau_DMP_res.res.acc

TC_DMP.t = TC_DMP_res.res.t
TC_DMP.p = TC_DMP_res.res.ref_pos
TC_DMP.v = TC_DMP_res.res.ref_vel
TC_DMP.a = TC_DMP_res.res.ref_acc

Opt_DMP.t = Opt_DMP_traj.traj.Time
Opt_DMP.p = Opt_DMP_traj.traj.Pos
Opt_DMP.v = Opt_DMP_traj.traj.Vel
Opt_DMP.a = Opt_DMP_traj.traj.Accel

%TODO: plot limits too!!


v_max = [2.1750; 2.1750; 2.1750; 2.1750; 2.6100; 2.6100; 2.6100];
a_max = [10; 7.5; 10; 10; 10; 10; 10];

figure('Name', 'Joint velocities')

for joint = 1:7
    if joint < 5
        subplot(2,4,joint)
    else
        subplot(2,4,joint+0.5)
    end

    if mod(joint,2) == 1
        tau_DMP.v(joint,:) = -tau_DMP.v(joint,:);
        TC_DMP.v(joint,:) = -TC_DMP.v(joint,:);
        Opt_DMP.v(joint,:) = -Opt_DMP.v(joint,:);
    end

    title(strcat('joint', num2str(joint), ' velocity'))
    hold on
    %plot_unconstrained = plot(unconstrained.t,unconstrained.v(joint,:),'m-', 'LineWidth', 1.5, 'DisplayName','unconstrained');
    plot_tau_DMP = plot(tau_DMP.t,tau_DMP.v(joint,:),'b-', 'LineWidth', 1.5, 'DisplayName','tau-DMP');
    plot_TC_DMP = plot(TC_DMP.t,TC_DMP.v(joint,:),'g-', 'LineWidth', 1.5, 'DisplayName','TC-DMP');
    plot_Opt_DMP = plot(Opt_DMP.t,Opt_DMP.v(joint,:),'k-', 'LineWidth', 1.5, 'DisplayName','Opt-DMP');

    limit_plot = plot([tau_DMP.t(1) tau_DMP.t(end)],[v_max(joint) v_max(joint)],'r:','LineWidth',1.5, 'DisplayName', 'Limit')
    plot([tau_DMP.t(1) tau_DMP.t(end)],-[v_max(joint) v_max(joint)],'r:','LineWidth',1.5)

    area([2,6],[v_max(joint),v_max(joint)], "FaceColor","y","FaceAlpha",0.3, "EdgeAlpha", 0)
    area([2,6],-[v_max(joint),v_max(joint)], "FaceColor","y","FaceAlpha",0.3, "EdgeAlpha", 0)

    legend([plot_tau_DMP, plot_TC_DMP, plot_Opt_DMP, limit_plot], 'Location','northwest')
end

figure('Name', 'Joint accelerations')
for joint = 1:7
    if joint < 5
        subplot(2,4,joint)
    else
        subplot(2,4,joint+0.5)
    end

    if mod(joint,2) == 1
        tau_DMP.a(joint,:) = -tau_DMP.a(joint,:);
        TC_DMP.a(joint,:) = -TC_DMP.a(joint,:);
        Opt_DMP.a(joint,:) = -Opt_DMP.a(joint,:);
    end

    title(strcat('joint', num2str(joint), ' acceleration'))
    hold on
    plot_tau_DMP = plot(tau_DMP.t,tau_DMP.a(joint,:),'k-', 'LineWidth', 1.5, 'DisplayName','tau-DMP');
    plot_TC_DMP = plot(TC_DMP.t,TC_DMP.a(joint,:),'g-', 'LineWidth', 1.5, 'DisplayName','TC-DMP');
    plot_Opt_DMP = plot(Opt_DMP.t,Opt_DMP.a(joint,:),'b-', 'LineWidth', 1.5, 'DisplayName','Opt-DMP');

    limit_plot = plot([tau_DMP.t(1) tau_DMP.t(end)],[a_max(joint) a_max(joint)],'r:','LineWidth',1.5, 'DisplayName', 'Limit')
    plot([tau_DMP.t(1) tau_DMP.t(end)],-[a_max(joint) a_max(joint)],'r:','LineWidth',1.5)

    legend([plot_tau_DMP, plot_TC_DMP, plot_Opt_DMP, limit_plot], 'Location','northwest')
end