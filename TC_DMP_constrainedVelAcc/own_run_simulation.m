function res = run_simulation(dmp,tmp_couple,params)

% Parameters
dt = params.dt;
T_max = params.T_max;

% Initialize
dmp = dmp.init();
t = 0;
k = 1;

while t < T_max
    % Store values
    res.t(k) = t;
    res.s(k) = dmp.s();
    res.pos(:,k) = dmp.ref_pos();
    res.vel(:,k) = dmp.ref_vel();
    res.acc(:,k) = dmp.ref_acc();
    res.tau(k) = dmp.tau;

    if norm(dmp.ref_pos()-dmp.g) < 1e-2
        break;
    end

    % Perform one time step
    k = k + 1;
    t = t + dt;
    
    % Step trajectory generator
    tau_dot = tmp_couple.tau_dot(dmp,dt);
    dmp = dmp.step(tau_dot,dt); 

end