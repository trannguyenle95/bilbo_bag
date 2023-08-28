function res = run_simulation(dmp,tmp_couple,params)

% Parameters
dt = params.dt;
T_max = params.T_max;
kp = params.kp;
kv = params.kv;
a_max = params.a_max;
v_max = params.v_max;
%s_v_max_online = params.s_v_max_online;
%v_max_online = params.v_max_online;

% Initialize
dmp = dmp.init();
t = 0;
k = 1;
sys.pos = dmp.ref_pos();
sys.vel = 0*sys.pos;
sys.acc = 0*sys.pos;

while t < T_max
    % Store values
    res.t(k) = t;
    res.s(k) = dmp.s();
    res.ref_pos(:,k) = dmp.ref_pos();
    res.ref_vel(:,k) = dmp.ref_vel();
    res.ref_acc(:,k) = dmp.ref_acc();
    res.sys_pos(:,k) = sys.pos;
    res.sys_vel(:,k) = sys.vel;
    res.sys_acc(:,k) = sys.acc;
    res.tau(k) = dmp.tau;

    %NOTE: might not need this? Can just let it run for T_max?
    % Final position convergence condition
    if norm(dmp.ref_pos()-dmp.g) < 1e-2
        break;
    end

    % Perform one time step
    k = k + 1;
    t = t + dt;
    
    %{
    % Add velocity limit online
    if dmp.s() < s_v_max_online
        v_max = v_max_online;
        tmp_couple.v_max= v_max_online;
    end
    %}

    % Step trajectory generator
    tau_dot = tmp_couple.tau_dot(dmp,dt);
    dmp = dmp.step(tau_dot,dt); 

    % Controller

    %NOTE: acceleration control with kp and kv ++ saturate a and v within
    %limits!! >> rewrite to just run reference DMP as is??
    
    %Possibility of using saturated controller for Franka? 
    %See rate limits here: https://frankaemika.github.io/docs/libfranka.html
    
    %{
    sys_acc_ctrl = dmp.ref_acc() + kv*(dmp.ref_vel()-sys.vel) + kp*(dmp.ref_pos()-sys.pos); %NOTE: DMP ref acc/vel/pos were updated in dmp.step() above!
    % Saturate acceleration
    sys_acc_ctrl = min(max(sys_acc_ctrl,-a_max),a_max);
    % Velocity integration
    vel_prev = sys.vel;
    sys.vel = sys.vel + sys_acc_ctrl*dt;
    % Saturate velocity
    sys.vel = min(max(sys.vel,-v_max),v_max);
    sys.acc = (sys.vel-vel_prev)/dt;
    % Position integration
    sys.pos = sys.pos + sys.vel*dt;
    %}
    
    %Direcly using DMP reference as system output
    sys.acc = dmp.ref_acc();
    sys.vel = dmp.ref_vel();
    sys.pos = dmp.ref_pos();
    %Now the system is not actually limited by vel and acc (no saturation)
    %so it looks like unconstrained DMP can follow motion (while breaking
    %vel and acc constraints) from plots too... >> this should be OK as we
    %can see from plots that constrained DMP now follows path while keeping
    %vel and acc wihtin bounds 

    %%NEXT: try on own reference trajs + 7 DOF (joints)
        %first some 2 DOF so plot works, just manage how to create struct
        %from csv
        %then all 7 joints + corresponding plots!


end