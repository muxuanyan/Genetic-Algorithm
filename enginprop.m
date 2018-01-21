function output = enginprop(L_1,L_2,L_3,y_m)
% Project #2 - Parametric Design of a MEMS Accelerometer
% Scott Moura
% SID 15905638
% ME 128, Prof. Lin
% Due Wed April 5, 2006
% Calculate Engineering Properties
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
% Define Inequality Constraint Limits
L_min = 20e-6; % Minimum beam length
a_r_min = 0.005*g; % Accerlation at minimum required resolution [m/s^2]
a_stress = 2000*g; % Acceleration at maximum required stress [m/s^2]
sigma_max = 1.6e9; % Maximum stress [Pa]
% Calculate design variable dependent parameters
x_m = L_b + 2*L_2 + 2*h_1 + 2*h_3; % Proof mass width [m]
A = x_m .* (y_m + 2*L_3 + 2*h_2); % Die Area [m^2]
zz = y_m + 2*L_3 + 2*h_2; % Die length [m]
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
sigma = norm([sigma_x sigma_y sigma_bend]);
% Find maximum survivable acceleration (assume only bending stress)
M_surv = sigma_max * I / (b_1/2);
s_surv = M_surv / (6*E*I .* (6*L_b*L_1.*L_3.^2 + 12*L_1.*L_2.*L_3.^2 + ...
2*L_b*L_2.*L_3.^2 + 6*L_1.*L_2.^2.*L_3 + L_b*L_2.^2.*L_3 + ...
4*L_b*L_1.*L_2.*L_3 + L_b*L_1.^2.*L_2) ./ LOCALden(L_1,L_2,L_3,L_b));
a_surv = (k ./ m) .* s_surv; % Maximum survivable acceleration
% Find displacements at 0.005g and 2000g
s_sm = m/k * 0.005*g;
s_lg = m/k * 2000*g;
% Find stress at displacement of 20 um
s_20 = 20e-6;
a_20 = (k ./ m) .* s_20;
M_20 = 6*s_20*E*I .* (6*L_b*L_1.*L_3.^2 + 12*L_1.*L_2.*L_3.^2 + ...
2*L_b*L_2.*L_3.^2 + 6*L_1.*L_2.^2.*L_3 + L_b*L_2.^2.*L_3 + ...
4*L_b*L_1.*L_2.*L_3 + L_b*L_1.^2.*L_2) ./ LOCALden(L_1,L_2,L_3,L_b);
F_x_20 = (k .* s_20) / 4;
F_y_20 = 18*s_20*E*I ./ L_2 .* (-L_b*L_1.^2.*L_2 + 2*L_b*L_1.*L_3.^2 - ...
2*L_b*L_1.^2.*L_3 + 6*L_1.*L_2.*L_3.^2 + L_b*L_2.*L_3.^2) ./ ...
LOCALden(L_1,L_2,L_3,L_b);
sigma_x_20 = F_x_20 / (h_1*b_1);
sigma_y_20 = F_y_20 / (h_1*b_1);
sigma_bend_20 = M_20 * (b_1/2) / I;
sigma_20 = norm([sigma_x_20 sigma_y_20 sigma_bend_20]);
% Find stress at 0.005g
M_sm = 6*s_sm*E*I .* (6*L_b*L_1.*L_3.^2 + 12*L_1.*L_2.*L_3.^2 + ...
2*L_b*L_2.*L_3.^2 + 6*L_1.*L_2.^2.*L_3 + L_b*L_2.^2.*L_3 + ...
4*L_b*L_1.*L_2.*L_3 + L_b*L_1.^2.*L_2) ./ LOCALden(L_1,L_2,L_3,L_b);
F_x_sm = (k .* s_sm) / 4;
F_y_sm = 18*s_sm*E*I ./ L_2 .* (-L_b*L_1.^2.*L_2 + 2*L_b*L_1.*L_3.^2 - ...
2*L_b*L_1.^2.*L_3 + 6*L_1.*L_2.*L_3.^2 + L_b*L_2.*L_3.^2) ./ ...
LOCALden(L_1,L_2,L_3,L_b);
sigma_x_sm = F_x_sm / (h_1*b_1);
sigma_y_sm = F_y_sm / (h_1*b_1);
sigma_bend_sm = M_sm * (b_1/2) / I;
sigma_sm = norm([sigma_x_sm sigma_y_sm sigma_bend_sm]);
% Display Engineering Properties
disp(['L_1: ' num2str(L_1*1e6) ' um'])
disp(['L_2: ' num2str(L_2*1e6) ' um'])
disp(['L_3: ' num2str(L_3*1e6) ' um'])
disp(['y_m: ' num2str(y_m*1e6) ' um'])
disp(['Die Area: ' num2str(A*1e12) ' um^2'])
disp(['Die Dimensions: ' num2str(x_m*1e6) ' um x ' num2str(zz*1e6) ' um'])
disp(['Proof Mass: ' num2str(m*1e9) ' ug'])
disp(['Spring Constant: ' num2str(k) ' N/m'])
disp(['Min Accel Resolution: ' num2str(a_res/g) 'g'])
disp(['Motion Resolution: ' num2str(s_res*1e6) ' um'])
disp(['Disp @ 0.005g: ' num2str(s_sm*1e6) ' um'])
disp(['Stress @ 0.005g: ' num2str(sigma_sm/1e9) ' GPa'])
disp(['Maximum Stress: ' num2str(sigma/1e9) ' GPA'])
disp(['Disp @ 2000g: ' num2str(s_lg*1e6) ' um'])
disp(['Max Survival Accel: ' num2str(a_surv/g) 'g'])
disp(['Disp @ ' num2str(a_surv/g) 'g: ' num2str(s_surv*1e6) ' um'])
disp(['Accel. @ 20 um disp: ' num2str(a_20/g) 'g'])
disp(['Stress @ 20 um disp: ' num2str(sigma_20/1e9) ' GPa'])