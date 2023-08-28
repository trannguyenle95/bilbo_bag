function D = ForwardKinematics(q)
    mdl_panda
    %info: https://www.petercorke.com/RTB/r9/html/SerialLink.html
    translation = panda.fkine(q);
    %https://mathworks.com/matlabcentral/fileexchange/83268-spatial-math-toolbox-for-matlab-peter-corke
    % >> see functions
    translation = translation.T;
    translation = translation(1:3,:,:); %trim away final row [0 0 0 1]

    pos = squeeze(translation(:,4,:))';
    quat = rotm2quat(translation(:,1:3,:)); %gives order [qw qx qy qz]
    quat = [quat(:,2) quat(:,3) quat(:,4) quat(:,1)]; %reordered to [qx qy qz qw]

    D = [pos quat];
end