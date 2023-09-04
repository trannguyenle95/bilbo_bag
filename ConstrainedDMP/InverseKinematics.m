function q = InverseKinematics(D)
    mdl_panda
    %info: https://www.petercorke.com/RTB/r9/html/SerialLink.html
    
    quats = [D(:,7) D(:,4) D(:,5) D(:,6)];
    T = quat2rotm(quats);
    T(:,4,:) = D(:,1:3)';
    T(4,:,:) = repmat([0 0 0 1], size(T, 3),1)';
    
    %For rotated gripper
    %q0 = [-0.06033718608193325, 0.1957925676774354, 0.1446464792309258, -2.1242161722067974, -0.12300981136857973, 2.3759525292393855, -0.6157846782301644];
    
    %Original gripper orientation
    %q0 = [-0.06153707569372121, 0.23268072162435294, -0.003733379824253959, -2.120620626949313, -0.07440938119840552, 2.374850448676014, 0.851590066155449];
    
    q0 = [0, 0.2837448589662732, 0, -2.0720574669683027, 0, 2.405712411822974, 0.7542077567525343]; %NOTE: position [0.59,0,0.20] gripper in original orientation


    q = panda.ikcon(T, q0);
    %q = panda.ikcon(T, q0, 'ObjectiveLimit', 1.0000e-04);
    %q = panda.ikine(T, q0, 'tol', 1e-4, 'alpha', 0.1);
end