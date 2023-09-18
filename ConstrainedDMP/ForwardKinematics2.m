function D = ForwardKinematics(q)
    franka = createFranka();
    translation = zeros(4,4,length(q));

    for i = 1:length(q)
        translation(:,:,i) = getTransform(franka,q(i,:)','body7');
    end

    translation = translation(1:3,:,:); %trim away final row [0 0 0 1]

    pos = squeeze(translation(:,4,:))';
    quat = rotm2quat(translation(:,1:3,:)); %gives order [qw qx qy qz]
    quat = [quat(:,2) quat(:,3) quat(:,4) quat(:,1)]; %reordered to [qx qy qz qw]

    D = [pos quat];
end