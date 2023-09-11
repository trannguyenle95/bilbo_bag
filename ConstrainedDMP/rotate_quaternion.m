function quat_rotated = rotate_quaternion(quats, axis, angle)
    %Use quaternion product here to rotate quaternions
    %https://en.wikipedia.org/wiki/Quaternion
    %first rotation q2 followed by rotation q1.

    qw2 = cos(angle/2);
    if axis == 'x'
        qx2 = sin(angle/2);
        qy2 = 0;
        qz2 = 0;
    elseif axis == 'y'
        qx2 = 0;
        qy2 = sin(angle/2);
        qz2 = 0;
    elseif axis == 'z'
        qx2 = 0;
        qy2 = 0;
        qz2 = sin(angle/2);
    else
        error("Invalid axis given for quaternion rotation")
    end
    
    qx1 = quats(:,1);
    qy1 = quats(:,2);
    qz1 = quats(:,3);
    qw1 = quats(:,4);
    
    quat_rotated(:,1) = qw1 .* qx2 + qx1 .* qw2 + qy1 .* qz2 - qz1 .* qy2;
    quat_rotated(:,2) = qw1 .* qy2 - qx1 .* qz2 + qy1 .* qw2 + qz1 .* qx2;
    quat_rotated(:,3) = qw1 .* qz2 + qx1 .* qy2 - qy1 .* qx2 + qz1 .* qw2;
    quat_rotated(:,4) = qw1 .* qw2 - qx1 .* qx2 - qy1 .* qy2 - qz1 .* qz2;
    
    quat_rotated(:,1:4) = quat_rotated(:,1:4) ./ sqrt(sum(quat_rotated(:,1:4).^2,2)); %Make sure quaternion is still unit

end