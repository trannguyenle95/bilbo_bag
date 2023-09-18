function jacobians = getJacobians(qs)
    jacobians = zeros(6,7,length(qs));
    franka = createFranka();
    for i = 1:length(qs)
        J = geometricJacobian(franka,qs(i,:)',"body7");
        jacobians(:,:,i) = J;
    end
end