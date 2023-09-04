mdl_panda
T = panda.fkine(q);

figure
hold on
plot(q(:,1), 'DisplayName','q1', 'LineWidth',2, 'Color', '#E1341E')
plot(q(:,2), 'DisplayName','q2', 'LineWidth',2, 'Color', '#13EC51')
plot(q(:,3), 'DisplayName','q3', 'LineWidth',2, 'Color', '#DB8A42')

plot(q(:,4), 'DisplayName','q4', 'LineWidth',2, 'Color',"#1ECBE1")
plot(q(:,5), 'DisplayName','q5', 'LineWidth',2, 'Color',"#EC13AE")
plot(q(:,6), 'DisplayName','q6', 'LineWidth',2, 'Color',"#4292DB")
plot(q(:,7), 'DisplayName','q7', 'LineWidth',2, 'Color',"#008F0B")

legend
hold off




M = T.double();

M(4,:,:) = []; %remove row with 0 0 1

t = squeeze(M(:,4, :))';
R = M(1:3, 1:3, :);

quats = rotm2quat(R); % qw qx qy qz


%get original data for comparison in plots

D = csvread("processed_sack_from_bag2.csv");


%plot(t)
figure
hold on
plot(t(:,1), 'DisplayName','x', 'LineWidth',5, 'Color', '#E1341E')
plot(t(:,2), 'DisplayName','y', 'LineWidth',5, 'Color', '#13EC51')
plot(t(:,3), 'DisplayName','z', 'LineWidth',5, 'Color', '#DB8A42')

plot(D(:,1), 'DisplayName','x_{in}', 'LineWidth',2, 'Color',"#1ECBE1")
plot(D(:,2), 'DisplayName','y_{in}', 'LineWidth',2, 'Color',"#EC13AE")
plot(D(:,3), 'DisplayName','z_{in}', 'LineWidth',2, 'Color',"#4292DB")

legend
hold off

figure
hold on
%Multiply by - to get same quaternion representation as in input
plot(-quats(:,1), 'DisplayName','qw', 'LineWidth',5,'Color', '#E1341E')
plot(-quats(:,2), 'DisplayName','qx', 'LineWidth',5,'Color', '#13EC51')
plot(-quats(:,3), 'DisplayName','qy', 'LineWidth',5,'Color', '#DB8A42')
plot(-quats(:,4), 'DisplayName','qz', 'LineWidth',5,'Color', '#6D927D')

plot(D(:,4), 'DisplayName','qx_{in}', 'LineWidth',2,'Color',"#EC13AE")
plot(D(:,5), 'DisplayName','qy_{in}', 'LineWidth',2,'Color',"#4292DB")
plot(D(:,6), 'DisplayName','qz_{in}', 'LineWidth',2,'Color', '#926D82')
plot(D(:,7), 'DisplayName','qw_{in}', 'LineWidth',2,'Color',"#1ECBE1")

legend
hold off
