function res = run_simulation(dmp,tmp_couple,params)

% Parameters
dt = params.dt;
T_max = params.T_max;

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
    
    % Step trajectory generator
    tau_dot = tmp_couple.tau_dot(dmp,dt);
    dmp = dmp.step(tau_dot,dt); 

    %Direcly using DMP reference as system output
    sys.acc = dmp.ref_acc();
    sys.vel = dmp.ref_vel();
    sys.pos = dmp.ref_pos();

end