% Project #2 - Parametric Design of a MEMS Accelerometer
% Scott Moura
% SID 15905638
% ME 128, Prof. Lin
% Due Wed April 5, 2006
% Genetic Algorithm
clear
% Define Equality Constraints
h_1 = 1.8e-6; % Beam Width [m]
h_2 = 1.8e-6; % Minimum Feature Length of Device
h_3 = 1.8e-6;
h_b = 1.8e-6;
b_1 = 1.8e-6; % Beam Thickness [m]
b_2 = 1.8e-6; % Minimum Feature Length of Device
b_3 = 1.8e-6;
b_m = 1.8e-6;
L_b = 150e-6; % Central Beam Width [m]
E = 160e9; % Elastic Modulus [Pa]
I = h_1*b_1^3/12; % Area Moment of Inertia [m^4]
rho = 2.33*100^3/1000; % Density [kg/m^3] (converted from gm/cm^3)
g = 9.8; % Gravitational Acceleration [m/s^2]
y_a = 50e-6; % Anchor Height [m]
x_a = L_b + 2*h_1; % Anchor Width [m]
A_max = 200000/(1e6)^2; % Maximum Die Area [m^2]
% Define Inequality Constraint Limits
L_min = 20e-6; % Minimum beam length
a_r_min = 0.005*g; % Accerlation at minimum required resolution [m/s^2]
a_stress = 2000*g; % Acceleration at maximum required stress [m/s^2]
sigma_max = 1.6e9; % Maximum stress [Pa]
% Define boundary limits
L_1_min = L_min;
L_1_max = 500e-6;
L_2_min = L_min;
L_2_max = 100e-6;
L_3_min = 100e-6;
L_3_max = 500e-6;
y_m_min = 100e-6;
y_m_max = 500e-6;
% Initialize best1 & best2
start_pops = 1000;
iters = 50;
best1 = zeros(iters,start_pops);
best2 = zeros(1,iters);
% Loop for 1000 different starting populations
for idx1 = 1:start_pops
% Define Set of Random Strings
d = 100; % Number of strings
clear Lambda
random = rand(d,4);
Lambda(:,1) = L_1_min + random(:,1) * (L_1_max - L_1_min); % L_1
Lambda(:,2) = L_2_min + random(:,2) * (L_2_max - L_2_min); % L_2
Lambda(:,3) = L_3_min + random(:,3) * (L_3_max - L_3_min); % L_3
Lambda(:,4) = y_m_min + random(:,4) * (y_m_max - y_m_min); % y_m
% Repeat Genetic Algorithm for 50 iterations
for idx2 = 1:iters
L_1 = Lambda(:,1);
L_2 = Lambda(:,2);
L_3 = Lambda(:,3);
y_m = Lambda(:,4);
% Calculate design variable dependent parameters
x_m = L_b + 2*L_2 + 2*h_1 + 2*h_3; % Proof mass width [m]
A = x_m .* (y_m + 2*L_3 + 2*h_2); % Die Area [m^2]
m = rho * x_m .* y_m * b_m; % Mass of system (ignore beams) [kg]
s_res = (1e-7)^2 ./ y_m; % Motion resolution [m]
% System Spring Constant [N/m]
k = 24*E*I * (6*L_1.*L_2.^2 + L_b*L_2.^2 + 4*L_b*L_1.*L_2 + 24*L_1.*L_2.*L_3 + ...
12*L_b*L_1.*L_3 + 4*L_b*L_2.*L_3) ./ LOCALden(L_1,L_2,L_3,L_b);
a_res = (k ./ m) .* s_res; % Minimum resolvable acceleration [m/s^2]
s_stress = (m ./ k) .* a_stress; % Displacement at maximum survival accel [m]
% Moment at max survival accel [N*m]
M = 6*s_stress*E*I .* (6*L_b*L_1.*L_3.^2 + 12*L_1.*L_2.*L_3.^2 + ...
2*L_b*L_2.*L_3.^2 + 6*L_1.*L_2.^2.*L_3 + L_b*L_2.^2.*L_3 + ...
4*L_b*L_1.*L_2.*L_3 + L_b*L_1.^2.*L_2) ./ LOCALden(L_1,L_2,L_3,L_b);
% Forces at max survival accel [N]
F_x = (k .* s_stress) / 4;
F_y = 18*s_stress*E*I ./ L_2 .* (-L_b*L_1.^2.*L_2 + 2*L_b*L_1.*L_3.^2 - ...
2*L_b*L_1.^2.*L_3 + 6*L_1.*L_2.*L_3.^2 + L_b*L_2.*L_3.^2) ./ ...
LOCALden(L_1,L_2,L_3,L_b);
% Stresses at max survival accel [N/m^2]
sigma_x = F_x / (h_1*b_1);
sigma_y = F_y / (h_1*b_1);
sigma_bend = M * (b_1/2) / I;
% Calculate the Fitness of the Objective Function
Pi = A*1e12;
% Check Constraints
for idx3 = 1:d
% Anchor cannot petrude through proof mass
if (L_3(idx3) < h_2 + y_a + L_1(idx3))
Pi(idx3) = inf;
end
% Check Die Area
if (A(idx3) > A_max)
Pi(idx3) = inf;
end
% Check minimum motion resolution
if (a_res(idx3) > a_r_min)
Pi(idx3) = inf;
end
% Check maximum stress
sigma_norm = norm([sigma_x(idx3) sigma_y(idx3) sigma_bend(idx3)]);
if (sigma_norm > sigma_max)
Pi(idx3) = inf;
end
% Check if beams crash into anchor
if (s_stress(idx3) > L_2(idx2))
%disp(s_stress(idx3))
Pi(idx3) = inf;
end
end
% Rank the Genetic Strings output by Pi and save the best
[sorted_Pi, order_Pi] = sort(Pi);
if (sorted_Pi(1) == inf)
%disp('INFINITY')
end
best1(idx2,idx1) = sorted_Pi(1);
% Keep the Ten Best Parents
nps = 10; % Number of Parents
for i = 1:nps
    parent(i,:) = Lambda(order_Pi(i),:);
end
% Mate the Top Ten Parents to Generate Ten Children
random2 = rand(nps,1);
for j = 1:2:nps-1
    child(j,:) = random2(j,1) * parent(j,:) + (1 - random2(j,1)) * parent(j+1,:);
    child(j+1,:) = random2(j+1,1) * parent(j,:) + (1 - random2(j+1,1)) * parent(j+1,:);
end
% Generate New Random Strings to Combine with Parents and Children
clear Lambda_new
random3 = rand(d-2*nps, 4);
Lambda_newrand(:,1) = L_1_min + random3(:,1) * (L_1_max - L_1_min); % new L_1
Lambda_newrand(:,2) = L_2_min + random3(:,2) * (L_2_max - L_2_min); % new L_2
Lambda_newrand(:,3) = L_3_min + random3(:,3) * (L_3_max - L_3_min); % new L_3
Lambda_newrand(:,4) = y_m_min + random3(:,4) * (y_m_max - y_m_min); % new y_m
Lambda_new = vertcat(parent, child, Lambda_newrand);
% Repeat Genetic Algorithm
Lambda = Lambda_new;
end
% Save the best lambda values
Lambda_final(idx1,:) = Lambda(1,:);
end
% Determine the runs that had the six best performing
best2 = best1(iters,:);
[sorted_best2, order_best2] = sort(best2);
% Output the design variable values for top five
for k = 1:5
    best_Lambda(k,:) = Lambda_final(order_best2(k),:);
    fprintf('Rank: %2.0f L_1 = %6.3e L_2 = %6.3e L_3 = %6.3e y_m = %6.3e Pi = %6.5f\n',...
    k,best_Lambda(k,1),best_Lambda(k,2),best_Lambda(k,3),best_Lambda(k,4),sorted_best2(k));
end
% Plot the minimum objective value as a function of generations
for idx4 = 1:100
    if (idx4 == 1)
        Pi_min(idx4) = best2(idx4);
    elseif (best2(idx4) <= Pi_min(idx4-1))
        Pi_min(idx4) = best2(idx4);
    else
        Pi_min(idx4) = Pi_min(idx4-1);
    end
end
plot(1:100,Pi_min(1:100))
title(['\bfMinimization of Objective Function: \Pi_{\it{min}} = ' num2str(sorted_best2(1))])
xlabel('Number of Starting Populations (First 100 Only)')
ylabel('Minimum Objective Value Found')