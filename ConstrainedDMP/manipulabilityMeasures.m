function [mu1, mu2, mu3, w, kappa] = manipulabilityMeasures(jacobians)
    mu1 = zeros(length(jacobians),1);
    mu2 = zeros(length(jacobians),1);
    mu3 = zeros(length(jacobians),1);
    w = zeros(length(jacobians),1);
    kappa = zeros(length(jacobians),1);

    for i = 1:length(jacobians)
        J = jacobians(:,:,i);
        
        eig_vals = eig(J*J');

        mu1(i) = sqrt(min(eig_vals)/max(eig_vals));
        mu2(i) = min(eig_vals)/max(eig_vals);
        %mu1(i) = sqrt(max(eig_vals)/min(eig_vals));
        %mu2(i) = max(eig_vals)/min(eig_vals);
        mu3(i) = sqrt(prod(eig_vals));
        w(i) = sqrt(det(J* J')); %alternative way to calculate mu3

        N1 = 1/size(J,2) * ones(size(J,2));
        J1 = sqrt(trace(J * N1 * J'));
        Jinv = pinv(J);
        N2 = 1/size(Jinv,2) * ones(size(Jinv,2));
        J2 = sqrt(trace(Jinv * N2 * Jinv'));
        kappa(i) = J1 * J2;
        
    end
end