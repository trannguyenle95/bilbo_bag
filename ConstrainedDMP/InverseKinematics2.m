function q = InverseKinematics2(D)
    q = zeros(length(D),7);

    franka = createFranka() %show(franka) <<< does it look reasonable?
    
    quats = [D(:,7) D(:,4) D(:,5) D(:,6)];
    T = quat2rotm(quats);
    T(:,4,:) = D(:,1:3)';
    T(4,:,:) = repmat([0 0 0 1], size(T, 3),1)';
    
    q0 = [0, 0.2837448589662732, 0, -2.0720574669683027, 0, 2.405712411822974, 0.7542077567525343]'; %NOTE: position [0.59,0,0.20] gripper in original orientation

    weights = [0.25 0.25 0.25 0.25 1 1];

    %IK only for first pose
    %Constraints defined in RigidBodyTree
    ik = inverseKinematics('RigidBodyTree',franka);

    q(1,:) = ik('body7',T(:,:,1),weights,q0);

    %[configSoln,solnInfo] = ik('panda_link7',pose1,weights,q0); %testing
    %imported panda

    %q(:,1) = configSoln;
    % NOTE - this worked before adding joint constraints but now returns "NaN"
    % in q when constraints have been added to joints in createFranka()
    
    for i = 2:length(q)
        q(i,:) = ik('body7',T(:,:,i),weights,q(i-1,:)');
    end
    



    %{
    gik = generalizedInverseKinematics('RigidBodyTree',franka, 'ConstraintInputs', {'jointbounds'});

    jointConst = constraintJointBounds(franka);
    jointConst.Bounds = [-2.7437, 2.7437;
                        -1.7628, 1.7628;
                        -2.8973, 2.8973;
                        -3.0421 -0.1518;
                        -2.8065, 2.8065;
                        0.5445, 3.7525;
                        -2.8973, 2.8973];

    configSol = gik(q0,jointConst);
    q= configSol
    %}

end