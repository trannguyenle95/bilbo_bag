function q = InverseKinematics2(D)
    q = zeros(length(D),7);
    
    %import franka
    franka = importrobot("frankaEmikaPanda.urdf")
    removeBody(franka, 'panda_rightfinger')
    removeBody(franka, 'panda_leftfinger')
    removeBody(franka, 'panda_hand')
    removeBody(franka, 'panda_link8')
    franka.DataFormat = 'column'; %so vector can be used instead of struct in inverse kinematics
    
    quats = [D(:,7) D(:,4) D(:,5) D(:,6)];
    T = quat2rotm(quats);
    T(:,4,:) = D(:,1:3)';
    T(4,:,:) = repmat([0 0 0 1], size(T, 3),1)';
    
    q0 = [0, 0.2837448589662732, 0, -2.0720574669683027, 0, 2.405712411822974, 0.7542077567525343]'; %NOTE: position [0.59,0,0.20] gripper in original orientation

    weights = [0.25 0.25 0.25 1 1 1];

    %IK only for first pose
    %Constraints defined in RigidBodyTree
    ik = inverseKinematics('RigidBodyTree',franka);

    q(1,:) = ik('panda_link7',T(:,:,1),weights,q0);
    
    for i = 2:length(q)
        q(i,:) = ik('panda_link7',T(:,:,i),weights,q(i-1,:)');
    end

end