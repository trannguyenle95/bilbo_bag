function [q, jacobians] = InverseKinematics2(D)
    q = zeros(length(D),7);
    
    jacobians = zeros(6,7,length(D));

    franka = createFranka();
    
    quats = [D(:,7) D(:,4) D(:,5) D(:,6)];
    T = quat2rotm(quats);
    T(:,4,:) = D(:,1:3)';
    T(4,:,:) = repmat([0 0 0 1], size(T, 3),1)';
    
    q0 = [0.1115568986356511, 0.4515083839684202, -0.06522623816080261, -2.4929104818055348, -0.11687397044566047, 2.9292743379890234, 0.8903726642525769]'; %for [0.45, 0, 0] and close to 45deg ee rotation


    weights = [0.25 0.25 0.25 1 1 1];  

    %Position constraints are defined in RigidBodyTree
    ik = inverseKinematics('RigidBodyTree',franka);
    q(1,:) = ik('body7',T(:,:,1),weights,q0);
    J = geometricJacobian(franka,q(1,:)',"body7");
    jacobians(:,:,1) = J;
    for i = 2:length(q)
        q(i,:) = ik('body7',T(:,:,i),weights,q(i-1,:)');
        J = geometricJacobian(franka,q(i,:)',"body7");
        jacobians(:,:,i) = J;
    end
    

end