clc;
clear;
close all;

%System Parameters
g = 9.81; l = 1; k = 1; m = 0.5;

%MPC Parameters
N = 20;             %prediction horizon steps
dt = 0.05;          %time step
Tsim = 5;           %total simulation time
Nsim = Tsim/dt;     %number of simulation steps

%Cost Weights
Q = diag([100, 1]);  %state cost (angle, velocity)
R = 0.1;             %input cost

%Constraints
u_max = 10;          %max torque
u_min = -10;

%Initial condition & goal
x = [0; 0];               %initial state
x_ref = [deg2rad(20); 0]; %desired upright position

%Storage for plotting
x_log = zeros(2, Nsim);
u_log = zeros(1, Nsim);

%Main Simulation Loop
for t = 1:Nsim
    %Solve NMPC optimization using fmincon
    u0 = zeros(N,1);  % initial guess
    lb = u_min*ones(N,1);
    ub = u_max*ones(N,1);

    cost_fun = @(u) nmpc_cost(u, x, x_ref, Q, R, dt, N, g, l, k, m);
    options = optimoptions('fmincon','Display','off');
    u_opt = fmincon(cost_fun, u0, [], [], [], [], lb, ub, [], options);

    %Apply first control input
    u = u_opt(1);
    x_dot = pendulum_dynamics(x, u, g, l, k, m);
    x = x + dt * x_dot;

    %Log data
    x_log(:,t) = x;
    u_log(t) = u;
end

%Plot Results
time = (0:Nsim-1)*dt;
figure;
subplot(3,1,1);
plot(time, rad2deg(x_log(1,:)),'LineWidth',1.5);
ylabel('Angle [deg]');
title('Nonlinear MPC - Pendulum');

subplot(3,1,2);
plot(time, rad2deg(x_log(2,:)),'LineWidth',1.5);
ylabel('Angular velocity [deg/s]');

subplot(3,1,3);
plot(time, u_log,'LineWidth',1.5);
ylabel('Torque u');
xlabel('Time [s]');

%Pendulum dynamics function
function dx = pendulum_dynamics(x, u, g, l, k, m)
    dx1 = x(2);
    dx2 = - (g/l)*sin(x(1)) - (k/m)*x(2) + (1/(m*l^2))*u;
    dx = [dx1; dx2];
end

%NMPC cost function
function J = nmpc_cost(u_seq, x0, x_ref, Q, R, dt, N, g, l, k, m)
    x = x0;
    J = 0;
    for i = 1:N
        u = u_seq(i);
        dx = pendulum_dynamics(x, u, g, l, k, m);
        x = x + dt * dx;
        e = x - x_ref;
        J = J + e'*Q*e + R*u^2;
    end
end
