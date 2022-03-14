%% Comparison of the RK4 solver as implemented, the fix, and ODE45 solution

% Assume zero control input
f_0 =  @(t,x)(-[ -2*x(2,:) ; 0.8*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:)] );
n = 2;
m = 1; % number of control inputs

x0 = [0.48, 0.25];

%% ************************ Original Discretization ***********************

deltaT = 0.01;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );        % k1 should be k3
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% ************************** Compare w/ ODE45 Solution ******************************

deltaT = 0.01;
ts = [0., 3.];
solODE45 = ode45(f_0, ts, x0);

T = []; X = [];
for t=0:deltaT:3.
    
end