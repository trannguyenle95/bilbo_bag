function D = preprocess(filename, displayplot, gripper_ori)

%D = readmatrix(filename)';
D = csvread(strcat('demos/',filename), 7); %offset 7 rows to skip header info

D = fillmissing(fillmissing(D, 'next'), 'previous'); %copy values if there are gaps, so quaternions should still be valid

if size(D, 2) == 9
    demoType = 'single'; %single hand tracked
else
    demoType = 'dual'; %two hands tracked
    x_dist = abs(D(:,7) - D(:,14)); %offset between both hands
end

%normalize around (0,0,0) add offset later
D(:,7) = D(:,7) - D(1,7); %x
D(:,8) = D(:,8) - D(1,8); %y
D(:,9) = D(:,9) - D(1,9); %y


%scale trajectory so height is within robot reach
max_h = max(D(:,8));
scale = 1;
if max_h > 0.6
    scale = 0.6 / max_h;
end
D(:,7:9) = D(:,7:9) * scale;


%assume rotation is mainly around x axis, if this is not the case the code
%must be edited
axis = 3;
D(:,7) = 0; %set linear movement around x axis to zero

%set rotation to be only around x axis
qx = D(:,axis); %copy qx
qw = sqrt(1-qx.^2) .* sign(D(:,6)); %make qw and qx add to 1, let qy and qz be 0

theta = 2*atan2(sqrt(qx.^2),qw); %calculate theta from new qx and qw

%offset rotation by -180 deg (-pi) because zero degree rotation means that
%the end-effector is pointing straight up (see visualization) and we want
%rotations to be with relation to the starting orientation where the
%gripper points straight downwards
qx = sin((theta - pi)/2);
qw = cos((theta - pi)/2); %-1 to change rotation direction

%insert new quaternion rotation, and also change order to q = [qw qx qy qz]
%notation from q = [qx qy qz qw] in OptiTrack
D(:, 3:6) = 0;
D(:, axis+1) = qx;
D(:, 3) = qw;
D = D(:,3:9); %skip index and time columns

%Change order of columns to match order (x,y,z,qx,qy,qz,qw) used by robot script
%and at the same times redefine axes to match change from OptiTrack frame
%to robot frame
reordered = D * 0; %create empty array
reordered(:,1) = D(:,5); %x
reordered(:,2) = D(:, 7) * -1; %y = old z flipped
reordered(:,3) = D(:, 6); %z = old y
reordered(:,4) = D(:, 2); %qx
reordered(:,5) = D(:, 4) * -1; %qy = old qz flipped
reordered(:,6) = D(:, 3); %qz = old qy
reordered(:,7) = D(:, 1); %qw
D = reordered;

%Now use quaternion rotation to rotate end-effector by -pi/2 around robot z-axis
%to match grip orientation.
%Use quaternion product here to rotate quaternions
%https://en.wikipedia.org/wiki/Quaternion
%first rotation q2 followed by rotation q1.
if gripper_ori == 1
    angle = +pi/4; %this gives sideways grip which needs higher distance between robot hands
    %q_rot = [cos(angle) 0 0 sin(angle)];
    %new_quat = quatmultiply(D(4:7),q_rot);
    qw2 = cos(angle/2);
    qx2 = 0;
    qy2 = 0;
    qz2 = sin(angle/2); %rotate around z axis
    
    qx1 = D(:,4);
    qy1 = D(:,5);
    qz1 = D(:,6);
    qw1 = D(:,7);
    
    D(:,4) = qw1 .* qx2 + qx1 .* qw2 + qy1 .* qz2 - qz1 .* qy2;
    D(:,5) = qw1 .* qy2 - qx1 .* qz2 + qy1 .* qw2 + qz1 .* qx2;
    D(:,6) = qw1 .* qz2 + qx1 .* qy2 - qy1 .* qx2 + qz1 .* qw2;
    D(:,7) = qw1 .* qw2 - qx1 .* qx2 - qy1 .* qy2 - qz1 .* qz2;
    
    D(:,4:7) = D(:,4:7) ./ sqrt(sum(D(:,4:7).^2,2)); %Make sure quaternion is still unit

else
    angle = -pi/4; %this gives "normal" grip which leaves more space between the grippers
    %q_rot = [cos(angle) 0 0 sin(angle)];
    %new_quat = quatmultiply(D(4:7),q_rot);
    qw2 = cos(angle/2);
    qx2 = 0;
    qy2 = 0;
    qz2 = sin(angle/2); %rotate around z axis
    
    qx1 = D(:,4);
    qy1 = D(:,5);
    qz1 = D(:,6);
    qw1 = D(:,7);
    
    D(:,4) = qw1 .* qx2 + qx1 .* qw2 + qy1 .* qz2 - qz1 .* qy2;
    D(:,5) = qw1 .* qy2 - qx1 .* qz2 + qy1 .* qw2 + qz1 .* qx2;
    D(:,6) = qw1 .* qz2 + qx1 .* qy2 - qy1 .* qx2 + qz1 .* qw2;
    D(:,7) = qw1 .* qw2 - qx1 .* qx2 - qy1 .* qy2 - qz1 .* qz2;
    
    D(:,4:7) = D(:,4:7) ./ sqrt(sum(D(:,4:7).^2,2)); %Make sure quaternion is still unit
end


% if gripper_ori == 2 %TEST
%     %Now use quaternion rotation to rotate end-effector by -pi/2 around robot z-axis
%     %to match grip orientation.
%     %Use quaternion product here to rotate quaternions
%     %https://en.wikipedia.org/wiki/Quaternion
%     %first rotation q2 followed by rotation q1.
%     angle = pi/4;
%     %q_rot = [cos(angle) 0 0 sin(angle)];
%     %new_quat = quatmultiply(D(4:7),q_rot);
%     qw2 = cos(angle/2);
%     qx2 = 0;
%     qy2 = 0;
%     qz2 = sin(angle/2); %rotate around z axis
% 
%     qx1 = D(:,4);
%     qy1 = D(:,5);
%     qz1 = D(:,6);
%     qw1 = D(:,7);
% 
%     D(:,4) = qw1 .* qx2 + qx1 .* qw2 + qy1 .* qz2 - qz1 .* qy2;
%     D(:,5) = qw1 .* qy2 - qx1 .* qz2 + qy1 .* qw2 + qz1 .* qx2;
%     D(:,6) = qw1 .* qz2 + qx1 .* qy2 - qy1 .* qx2 + qz1 .* qw2;
%     D(:,7) = qw1 .* qw2 - qx1 .* qx2 - qy1 .* qy2 - qz1 .* qz2;
% 
%     D(:,4:7) = D(:,4:7) ./ sqrt(sum(D(:,4:7).^2,2)); %Make sure quaternion is still unit
% end


if strcmp(demoType,'dual')
    %add offsets
    %robots are 1.5m apart, so for total distance "x_dist" the distance from the
    %middle point (0.75m) to either robot should be "x_dist"/2.
    if gripper_ori == 1
        %remember wide side ~20cm
        %Gripper with wide side towards eachother = smaller max distance 
        %min distance between grippers is (0.75-0.575)*2 - 0.20 = 0.15 between outer parts of grippers.
        %Set maximum distance between grippers to (0.75 - 0.50)*2 = 0.50.
        D(:,1) = max(min(0.75 - (x_dist/2), 0.575),0.50); %x
    else
        %short side of grippers ~10cm
        %Gripper with short side towards eachother = larger max distance 
        %min distance between grippers is (0.75-0.625)*2 - 0.10 = 0.15 between outer parts of grippers.
        %Set maximum distance between grippers to (0.75 - 0.50)*2 = 0.50.
        %D(:,1) = max(min(0.75 - (x_dist/2), 0.625),0.50); %x - made DMP diverge from demo in x direction
        %D(:,1) = max(min(0.75 - (x_dist/2), 0.60),0.50); %x - made DMP diverge from demo in x direction for other demo
        D(:,1) = max(min(0.75 - (x_dist/2), 0.575),0.50); %x
    end
else
    if gripper_ori == 1
        D(:,1) = 0.575; %x
    else
        %D(:,1) = 0.625; %x - made DMP diverge from demo in x direction
        %D(:,1) = 0.60; %x - made DMP diverge from demo in x direction
        D(:,1) = 0.575; %x
    end
end

D(:,2); %y
D(:,3) = D(:,3) + 0.20; %z

%Flip y direction
%D(:,2) = -D(:,2); %position
%D(:,5) = -D(:,5); %rotation - is zero unless ee is rotated

%%Flip rotation around x!
%D(:,4) = -D(:,4); %flip rotation direction

%plot
if displayplot == true
    figure()
    hold on
    plot(D(:,1),'r-','LineWidth',1.5, 'DisplayName','x');
    plot(D(:,2),'g-','LineWidth',1.5, 'DisplayName','y');
    plot(D(:,3),'b-','LineWidth',1.5, 'DisplayName','z');
    legend
    hold off
    
    figure()
    hold on
    plot(D(:,4),'r-','LineWidth',1.5, 'DisplayName','qx');
    plot(D(:,5),'g-','LineWidth',1.5, 'DisplayName','qy');
    plot(D(:,6),'b-','LineWidth',1.5, 'DisplayName','qz');
    plot(D(:,7),'k-','LineWidth',1.5, 'DisplayName','qw');
    legend
    hold off
end

end