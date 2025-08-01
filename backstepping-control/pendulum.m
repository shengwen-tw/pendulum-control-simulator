clc; clear; close all;

%System Parameters
g = 9.81;   %gravity
l = 1;      %length
k = 1;      %damping coefficient
m = 0.5;    %mass

%Backstepping gains
k1 = 5;     %for e1
k2 = 5;     %for e2

%Simulation Setup
dt = 0.001;
T = 10;
N = T / dt;

x1 = deg2rad(0);      %initial angle
x1_des = deg2rad(20); %desired angle
x2 = 0;               %initial angular velocity

x1_arr = zeros(1,N);
x2_arr = zeros(1,N);
u_arr  = zeros(1,N);
t_arr  = zeros(1,N);

for i = 1:N
    %Errors
    e1 = x1 - x1_des;
    x2_des = -k1 * e1;
    e2 = x2 - x2_des;
    
    %Control Law (Backstepping)
    u = m*l^2 * ( -k1*x2 - k2*e2 + (g/l)*sin(x1) + (k/m)*x2 );
    
    %Dynamics
    x1_dot = x2;
    x2_dot = (-g/l)*sin(x1) - (k/m)*x2 + (1/(m*l^2)) * u;
    
    %Integration
    x1 = x1 + x1_dot * dt;
    x2 = x2 + x2_dot * dt;

    %Logging
    x1_arr(i) = rad2deg(x1);
    x2_arr(i) = rad2deg(x2);
    u_arr(i) = u;
    t_arr(i) = i * dt;
end

%Plot
figure;
subplot(3,1,1);
plot(t_arr, x1_arr, 'LineWidth', 1.5);
ylabel('Angle [deg]');
title('Backstepping Control of Inverted Pendulum');
grid on;

subplot(3,1,2);
plot(t_arr, x2_arr, 'LineWidth', 1.5);
ylabel('Angular Velocity [deg/s]');
grid on;

subplot(3,1,3);
plot(t_arr, u_arr, 'LineWidth', 1.5);
ylabel('Control Input');
xlabel('Time [s]');
grid on;

pause;
