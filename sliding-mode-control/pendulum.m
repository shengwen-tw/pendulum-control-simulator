clc; clear; close all;

%Initial State
x1 = deg2rad(0);  % Initial angle (rad)
x2 = 0;           % Initial angular velocity (rad/s)
u = 0;            % Initial control input

%System Parameters
g = 9.8;          % gravity (m/s^2)
l = 1;            % rod length (m)
k = 1;            % air drag coefficient
m = 0.5;          % mass (kg)

%Sliding Mode Control Parameters
c = 5;            % sliding surface slope
eta = 10;         % switching gain (larger = faster but more chattering)

%Desired Setpoint
x1_d = deg2rad(100);   % desired angle
x2_d = 0;              % desired angular velocity
x_d = [x1_d; x2_d];

%Simulation Time
ITERATION_TIMES = 10000;
dt = 0.001;  %1ms step, total 10s

%Data Logging
time_arr = zeros(1, ITERATION_TIMES);
x1_arr = zeros(1, ITERATION_TIMES);
x2_arr = zeros(1, ITERATION_TIMES);
u_arr  = zeros(1, ITERATION_TIMES);
s_arr  = zeros(1, ITERATION_TIMES);

%Simulation Loop
tic();
for i = 1:ITERATION_TIMES
    %System dynamics update
    x1_dot = x2;
    x2_dot = (-g/l)*sin(x1) - (k/m)*x2 + (1/(m*l*l)) * u;
    x1 = x1 + x1_dot * dt;
    x2 = x2 + x2_dot * dt;
    x = [x1; x2];

    %Sliding Mode Control
    e = x - x_d;
    s = c * e(1) + e(2); %sliding surface
    u = (m*l^2) * ( -c*x2 + (g/l)*sin(x1) + (k/m)*x2 - eta * sign(s) );

    %Logging
    x1_arr(i) = rad2deg(x1);
    x2_arr(i) = rad2deg(x2);
    time_arr(i) = i * dt;
    u_arr(i) = u;
    s_arr(i) = s;
end
toc();

%Plotting
figure;
subplot(3,1,1);
plot(time_arr, x1_arr, 'LineWidth', 1.5);
title('Pendulum Angle (Sliding Mode Control)');
ylabel('Angle [deg]');
grid on;

subplot(3,1,2);
plot(time_arr, x2_arr, 'LineWidth', 1.5);
ylabel('Angular Velocity [deg/s]');
grid on;

subplot(3,1,3);
plot(time_arr, u_arr, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Control Input u');
grid on;

%Plot Sliding Surface
figure;
plot(time_arr, s_arr, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Sliding Surface s');
title('Sliding Surface Convergence');
grid on;

pause;
