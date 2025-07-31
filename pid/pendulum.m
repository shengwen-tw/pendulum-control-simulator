clear;
clc;

%System parameters
g = 9.8;
l = 1;
k = 1;
m = 0.5;

%PID parameters
Kp = 30;
Ki = 10;
Kd = 5;

%Reference step
theta_ref_deg = 100; % step input
theta_ref_rad = deg2rad(theta_ref_deg);

%s-domain analysis
s = tf('s');
G = 2 / (s^2 + 2*s + 9.8);
C = Kp + Ki/s + Kd*s;
T = feedback(C*G, 1);

%Step response + performance metrics
figure(1);
step(T);
title('Step Response of Linearized System (s-domain)');
ylabel('Angle [rad]');
xlabel('Time [s]');
grid on;

%Step response info
info = stepinfo(T);
closed_loop_poles = pole(T);
is_stable = all(real(closed_loop_poles) < 0);

%Print performance with explanations
fprintf('\n[S-domain linear system performance]\n');

fprintf('Overshoot: %.2f %%\n', info.Overshoot);
fprintf('  → The maximum percentage the response exceeds the reference value.\n');
fprintf('     Formula: (Peak - Target) / Target × 100%%\n\n');

fprintf('Rise Time: %.4f s\n', info.RiseTime);
fprintf('  → Time it takes for the response to rise from 10%% to 90%% of the target.\n\n');

fprintf('Settling Time: %.4f s\n', info.SettlingTime);
fprintf('  → Time for the response to stay within ±2%% of the target without leaving again.\n\n');

fprintf('All poles in LHP (stable)? %s\n', mat2str(is_stable));
fprintf('  → System is stable if all poles have negative real parts (left half-plane).\n\n');

fprintf('Closed-loop poles:\n');
disp(closed_loop_poles);

%Time-domain simulation (nonlinear)
x1 = 0;
x2 = 0;
u = 0;
dt = 0.001;
ITERATION_TIMES = 10000;

time_arr = zeros(1, ITERATION_TIMES);
x1_arr = zeros(1, ITERATION_TIMES);
x2_arr = zeros(1, ITERATION_TIMES);

integral_error = 0;
prev_error = 0;

for i = 1:ITERATION_TIMES
    t = (i - 1) * dt;

    % nonlinear pendulum dynamics
    x1_dot = x2;
    x2_dot = (-g/l)*sin(x1) - (k/m)*x2 + (1/(m*l^2))*u;
    x1 = x1 + x1_dot * dt;
    x2 = x2 + x2_dot * dt;

    % PID controller
    error = theta_ref_rad - x1;
    integral_error = integral_error + error * dt;
    derivative_error = (error - prev_error) / dt;
    prev_error = error;

    u = Kp * error + Ki * integral_error + Kd * derivative_error;

    % log
    x1_arr(i) = x1;
    x2_arr(i) = x2;
    time_arr(i) = t;
end

%Plot nonlinear result
figure(2);
subplot(2, 1, 1);
plot(time_arr, rad2deg(x1_arr));
title('Nonlinear Pendulum Simulation (time domain)');
xlabel('Time [s]');
ylabel('Angle [deg]');
grid on;

subplot(2, 1, 2);
plot(time_arr, rad2deg(x2_arr));
xlabel('Time [s]');
ylabel('Angular Velocity [deg/s]');
grid on;

%Pole-Zero Map
figure(3); clf;
pzmap(T);
title('Pole-Zero Map of Closed-Loop System');
grid on;
