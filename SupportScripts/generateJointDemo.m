%this is done in each DMP version so keep it in a shared script

if strcmp('A', bag)
    bagwidth = 0.38
elseif strcmp('B', bag)
    bagwidth = 0.40
elseif strcmp('C', bag)
    bagwidth = 0.42
elseif strcmp('D', bag)
    bagwidth = 0.49
elseif strcmp('E', bag)
    bagwidth = 0.53
end

D = preprocess(filename, false, 0.00, 0.05, 0.00, 1, 'ori1', bagwidth);
Dsmooth = smoothdata(D, 1, "gaussian", 35); %smooth demo before calculating IK
Dsmooth(:,4:7) = Dsmooth(:,4:7) ./ sqrt(sum(Dsmooth(:,4:7).^2,2)); %Make sure quaternion still has unit norm
[q, jacobians] = InverseKinematics(Dsmooth);
demo_traj = generateDemo(q', 1/120);