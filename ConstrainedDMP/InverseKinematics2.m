function [q, jacobians] = InverseKinematics2(D)
    q = zeros(length(D),7);
    
    jacobians = zeros(6,7,length(D));

    franka = createFranka();
    
    quats = [D(:,7) D(:,4) D(:,5) D(:,6)];
    T = quat2rotm(quats);
    T(:,4,:) = D(:,1:3)';
    T(4,:,:) = repmat([0 0 0 1], size(T, 3),1)';
    
    %q0 = [0, 0.2837448589662732, 0, -2.0720574669683027, 0, 2.405712411822974, 0.7542077567525343]'; %NOTE: position [0.59,0,0.20] gripper in original orientation
    q0 = [0.1115568986356511, 0.4515083839684202, -0.06522623816080261, -2.4929104818055348, -0.11687397044566047, 2.9292743379890234, 0.8903726642525769]'; %for [0.45, 0, 0] and close to 45deg ee rotation
    %q0 = [0.11424330942149746, 0.6734107206244216, -0.08293315969712552, -1.988764453510407, -0.1179131529739628, 2.655966289573245, 0.7821851946132062]'; %for [0.60, 0, 0] and close to 45deg ee rotation
    %q0 = [0.16169800343829296, 0.2535896401530817, -0.10843374721545593, -2.938105589615671, -0.12041141999328904, 3.2314470050070017, 0.7651126211285408]'; %for [0.30, 0, 0] and close to 45deg ee rotation
    %q0 = [0.04705703931386988, 0.9539009694132887, -0.036918611843170177, -2.210548720291037, -0.11555665776050751, 3.180519715680016, 0.8974417739769928]'; %close to [0.45, 0, -0.15]
    %q0 = [0.05605581243519197, 0.060367270439538136, -0.04229283738815993, -2.5147607649070762, -0.11617064242561657, 2.6120584942616585, 0.892399529736696]'; %close to [0.45, 0, 0.15]


    weights = [0.25 0.25 0.25 1 1 1];   

    %IK only for first pose
    %Constraints defined in RigidBodyTree
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