function q = InverseKinematics(D)
    mdl_panda
    %info: https://www.petercorke.com/RTB/r9/html/SerialLink.html
    
    quats = [D(:,7) D(:,4) D(:,5) D(:,6)];
    T = quat2rotm(quats);
    T(:,4,:) = D(:,1:3)';
    T(4,:,:) = repmat([0 0 0 1], size(T, 3),1)';
    
    q0 = [-0.06033718608193325, 0.1957925676774354, 0.1446464792309258, -2.1242161722067974, -0.12300981136857973, 2.3759525292393855, -0.6157846782301644];
    q = panda.ikcon(T, q0);
    %q = panda.ikcon(T, q0, 'ObjectiveLimit', 1.0000e-04);
    %q = panda.ikine(T, q0, 'tol', 1e-4, 'alpha', 0.1);
end