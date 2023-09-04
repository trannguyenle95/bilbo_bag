q_vel = diff(q) / (1/1000);

q_acc = diff(q_vel) / (1/1000);

figure
hold on
plot(q_vel(:,1), 'DisplayName','v_{q1}')
plot(q_vel(:,2), 'DisplayName','v_{q2}')
plot(q_vel(:,3), 'DisplayName','v_{q3}')
plot(q_vel(:,4), 'DisplayName','v_{q4}')
plot(q_vel(:,5), 'DisplayName','v_{q5}')
plot(q_vel(:,6), 'DisplayName','v_{q6}')
plot(q_vel(:,7), 'DisplayName','v_{q7}')
legend
hold off

figure
hold on
plot(q_acc(:,1), 'DisplayName','a_{q1}')
plot(q_acc(:,2), 'DisplayName','a_{q2}')
plot(q_acc(:,3), 'DisplayName','a_{q3}')
plot(q_acc(:,4), 'DisplayName','a_{q4}')
plot(q_acc(:,5), 'DisplayName','a_{q5}')
plot(q_acc(:,6), 'DisplayName','a_{q6}')
plot(q_acc(:,7), 'DisplayName','a_{q7}')
legend
hold off