%% Comparison of the RK4 solver as implemented, the fix, and ODE45 solution
clear; close all; clc;

u = 0.5;
% Assume zero control input
f_u =  @(t,x,u)(-[ -2*x(2,:) ; 0.8*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:) + u] );
n = 2;
m = 1; % number of control inputs

X0 = [0.48; 0.25];

%% ************************ Original Discretization ***********************

deltaT = 0.01;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );        % k1 should be k3
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% ************************ Fixed Discretization ***********************

deltaT = 0.01;
%Runge-Kutta 4
k4_fix = @(t,x,u) ( f_u(t,x + k3(t,x,u)*deltaT,u) );        % k1 should be k3
f_ud_fix = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4_fix(t,x,u)  )   );

%% ************************** Compare w/ ODE45 Solution ******************************

deltaT = 0.01;
ts = [0., 10.];
solODE45 = ode45(@(t,x) f_u(t,x,u), ts, X0);

T = []; X = []; X_fix = [];
Xk = X0;
Xk_fix = X0;
for t=0:deltaT:3.
    T = [T t];
    X = [X Xk];
    X_fix = [X_fix Xk_fix];
    Xk = f_ud(deltaT, Xk, u);
    Xk_fix = f_ud_fix(deltaT, Xk_fix, u);
end

%% Plot the results
figure; plot(T,X(1,:), 'b--', 'LineWidth', 2); hold on;
plot(T,X_fix(1,:), 'r-', 'LineWidth', 1);
plot(solODE45.x, solODE45.y(1,:), 'k.', 'LineWidth',2);

figure; plot(T,X(2,:), 'b--', 'LineWidth', 2); hold on;
plot(T,X_fix(2,:), 'r-', 'LineWidth', 1);
plot(solODE45.x, solODE45.y(2,:), 'k.', 'LineWidth',2);