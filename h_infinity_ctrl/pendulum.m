x1 = deg2rad(0);  %angle
x2 = 0;           %angular velocity
u = 0;            %torque, control input
d = 0;            %torque, external disturbance
g = 9.8;          %gravitational acceleration
l = 1;            %length of the rod [m]
k = 1;            %air drag coefficient
m = 0.5;          %mass [kg]

%controller setpoint
x1_d = deg2rad(100); %desired angle
x2_d = deg2rad(0); %desired angular velocity
x_d = [x1_d; x2_d];

%simulation time: ITERATION_TIMES * dt = ~10s
ITERATION_TIMES = 10000;
dt = 0.001;

%array for plotting
time_arr = zeros(1, ITERATION_TIMES);
x1_arr = zeros(1, ITERATION_TIMES);
x2_arr = zeros(1, ITERATION_TIMES);
d_arr = zeros(1, ITERATION_TIMES);
gamma_arr = zeros(1, ITERATION_TIMES);

tic();
for i = 1: ITERATION_TIMES
    disp(i);
    
    %=============================%
    % generate random disturbance %
    %=============================%
    sigma_d = 360; %distribution of the torque disturbance
    tau_c = 3.2;   %correlation time of the disturbance ODE
    %
    inv_cor_time = -1 / tau_c;
    A_d = (-1 / tau_c) .* eye(1, 1);
    random_noise = sigma_d * randn(1, 1);
    d_dot = (A_d * d) + random_noise;
    d = d + d_dot * dt;
    
    %========================%
    % update system dynamics %
    %========================%
    x1_dot = x2;
    x2_dot = (-g/l)*sin(x1) - (k/m)*x2 + (1/(m*l*l))*u + (1/(m*l*l))*d;
    x1 = x1 + x1_dot * dt;
    x2 = x2 + x2_dot * dt;
    x = [x1; x2];

    %====================%
    % H-infinity control %
    %====================%
    
    %state transition matrix
    a11 = 0;
    a12 = 1;
    a21 = (-g/l)*cos(x1);
    a22 = -k/m;
    A = [a11 a12;
         a21 a22];
    
    %===========================%
    % x_dot = A*x + B1*w + B2*u %
    %===========================%
    
    %disturbance input matrix
    B1 = [0;
          1/(m*l*l)];
      
    %control input matrix
    B2 = [0;
          1/(m*l*l)];
    
    %==================%
    % z = C1*x + D12*u %
    %==================%
    
    %penalty error measurement (weighting) matrix
    C1 = [600 0   %weight of minimizing the angle error
          0   5   %weight of minimizing the angular velocity error
          0   0];
      
    %penalty error control input (weighting) matrix
    w_u = 0.1;        %weight of the control energy
    D12 = [0;
           0;
           w_u/(m*l*l)];
       
    %==================%
    % y = C2*x + D21*w %
    %==================%
    
    %state-feedback H-infinity control problem, i.e, C2 = I, D21 = 0
    C2 = eye(2);
    D21 = 0;
    
    %====================%
    % H-infinity control %
    %====================%
    [gamma, X] = hinf_syn(A, B1, B2, C1, 0);
    
    %calculate the feedback control
    C0_hat = -B2.' * X;
    u = C0_hat * [x - x_d];
    
    %sys=ss(A + B2*C0_hat, B1, C1 + D12*C0_hat, 0);
    %sigma(sys, ss(gamma));
    %pause;
    
    %record datas for plotting
    x1_arr(i) = rad2deg(x1);
    x2_arr(i) = rad2deg(x2);
    d_arr(i) = d;
    gamma_arr(i) = gamma;
    time_arr(i) = i * dt;
end
toc();

subplot (2, 1, 1);
plot(time_arr, x1_arr);
title('Pendulum H-infinity Control');
xlabel('time [s]');
ylabel('angle [deg]');
subplot (2, 1, 2);
plot(time_arr, x2_arr);
xlabel('time [s]');
ylabel('angular velocity [deg/s]');

figure;
plot(time_arr, gamma_arr);
title('\gamma (||G_{\infty}||)');
xlabel('time [s]');
ylabel('\gamma (||G_{\infty}||)');

figure;
plot(time_arr, d_arr);
title('Torque disturbance');
xlabel('time [s]');
ylabel('\tau_d');

disp("Press any key to leave");
pause;
close all;