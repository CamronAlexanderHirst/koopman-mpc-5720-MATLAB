clear all
close all
addpath('./Resources')
rng(2141444)

% This script performs Koopman linearization of a simplified Dubin's AC
% model, where the position is invariant. All user-adjustable parameters
% are located at the top of the script with descriptions.

%%% Challenges %%%
% The current linearization uses 1000 datapoints and 4 
% basis fucntions. The challenge is to use a combination of inverse
% quadratic and thin plate basis functions to get a "good" linearization
% with as little data as possible. 

% 1. Can you adjust the mixture of basis functions to improve the
% prediction quality (Lines 52 and 57)? You may need 
% more data points to accomplish this (Lines 42-43). 

% 2. Increase the predicition horizon (Tmax Line 68) - how far out is the 
% Koopman approximation good? Can you get arbitrarily close prediction? 
% Does it work for other initial conditions (x0 Line 75)?

% 3. Change the control signal to the constant signal (u_dt Line 71), and 
% make sure the regression is still good and make changes if necessary.

% 4. Try out the random control signal -- can you find a choice of basis,
% data, etc. to get good predictions?

% 5. Are there any basis functions that might make sense besides these
% radial basis functions for this system? Add functions to the basis 
% function vector to try them out (Line 63).

% 6. Try defining your own input signals or dynamics (Line 84) - you may 
% have to change the training data range to get it to work the way you 
% want (Lines 105 and 108)!

%%%%%%%%%%%%%%%%%%%%%%%%%%% Data collection parameters

Nsim = 10;                    % # of simulations
Ntraj = 100;                   % # of steps in each simulation
heading_sample_range = pi;     % +/- this range

%%%%%%%%%%%%%%%%%%%%%%%%%%% Basis functions

% Basis Function Definitions
n = 3;                                % Dimension of system (do not change)

% RBF 1 Centers
Nrbf = 2;                             % # of basis functions
cent = rand(n,Nrbf)*2 - 1;            % centers of each function
rbf_type = 'invquad'; 

% RBF 2 Centers
Nrbf2 = 2;                            % # of basis functions
cent2 = rand(n,Nrbf2)*2 - 1;          % centers of each function
rbf_type_2 = 'thinplate';             % type of function - one of 'thinplate', 'gauss', 'invquad', 'polyharmonic'

% Lifting mapping - RBFs + the state itself
extra_param = 1;
liftFun = @(xx)( [
                 xx; 
                 rbf(xx(3,:),cent,rbf_type,extra_param); 
                 rbf(xx(3,:), cent2,rbf_type_2);
                 ]);
Nlift = (Nrbf + Nrbf2) + n;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation Parameters
% Simulation Results
Tmax = 1;                   % Max Timestep

% Control signal
u_dt = @(i)((-1).^(round(i/250))); % Stepping signal
% u_dt = @(i)(0.5);                % Constant signal
% u_dt = @(i) 2*rand(1)-1;

x0 = [0.; 0.; pi/4];       % Initial Condition (North / East / Heading)  
        
%% *************************** Dynamics ***********************************
v = 1.;
% Assume zero control input
f_u =  @(t,x,u)([ v*cos(x(3,:)) ; v*sin(x(3,:)) ; u ] );

n = 3;
m = 1; % number of control inputs

%% ************************** Discretization ******************************

deltaT = 0.01;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k3(t,x,u)*deltaT,u) ); 
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% ************************** Collect data ********************************
tic
disp('Starting data collection')

% Random forcing
Ubig = 2*rand([Nsim Ntraj]) - 1;

% Random initial heading
Xcurrent = [zeros(2,Ntraj); mod(rand(1,Ntraj)*heading_sample_range*2 - heading_sample_range, 2*pi)];


X = []; Y = []; U = [];
for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Ubig(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Xcurrent = Xnext;
end
fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ******************************* Lift ***********************************

disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);
fprintf('Lifting DONE, time = %1.2f s \n', toc);

%% ********************** Build predictor *********************************

disp('Starting REGRESSION')
tic
W = [Ylift ; X];
V = [Xlift; U];
VVt = V*V';
WVt = W*V';
M = WVt * pinv(VVt); % Matrix [A B; C 0]
Alift = M(1:Nlift,1:Nlift);
Blift = M(1:Nlift,Nlift+1:end);
Clift = M(Nlift+1:end,1:Nlift);

fprintf('Regression done, time = %1.2f s \n', toc);

%% *********************** Predictor comparison ***************************

Nsim = Tmax/deltaT;

% Initial condition
x_true = x0;

% Lifted initial condition
xlift = liftFun(x0);

% Local linearization predictor at x0
x = sym('x',[3;1]); u = sym('u',[1;1]);

Ac_x0 = double(subs(jacobian(f_u(0,x,u),x),[x;u],[x0;0]));
Bc_x0 = double(subs(jacobian(f_u(0,x,u),u),[x;u],[x0;0]));
c_x0 = double(subs(f_u(0,x,u),[x;u],[x0;0])) - Ac_x0*x0 - Bc_x0*0;
ABc = expm([Ac_x0 [Bc_x0 c_x0] ; zeros(2,5)]*deltaT); % discretize
Ad_x0 = ABc(1:3,1:3); Bd_x0 = ABc(1:3,3); cd_x0 = ABc(1:3,4);
X_loc_x0 = x0;

% % Local linearization predictor at 0
x = sym('x',[3;1]); u = sym('u',[1;1]);
Ac_0 = double(subs(jacobian(f_u(0,x,u),x),[x;u],[0;0;0;0]));
Bc_0 = double(subs(jacobian(f_u(0,x,u),u),[x;u],[0;0;0;0]));
c_0 = double(subs(f_u(0,x,u),[x;u],[0;0;0;0])) - Ac_0*[0;0;0] - Bc_0*0;
ABc = expm([Ac_0 [Bc_0 c_0] ; zeros(2,5)]*deltaT); % discretize
Ad_0 = ABc(1:3,1:3); Bd_0 = ABc(1:3,3); cd_0 = ABc(1:3,4); 
X_loc_0 = x0;

X_koop = [x0];
X_loc_new = [x0];
% Simulate
i_relinearize = 30;
for i = 0:Nsim-1
    x_tmp = [0; 0; X_koop(3,end)];
    x_tmp_lift = liftFun(x_tmp);
    delta_x_lift = Alift*x_tmp_lift + Blift*u_dt(i);
    X_koop = [X_koop [X_koop(1:2,end)+delta_x_lift(1:2); mod(delta_x_lift(3), 2*pi)]];
    
    % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
    
%     % Local linearization predictor at x0
    X_loc_x0 = [X_loc_x0, Ad_x0*X_loc_x0(:,end) + Bd_x0*u_dt(i) + cd_x0];

%     % Local linearization predictor at 0
    X_loc_0 = [X_loc_0, Ad_0*X_loc_0(:,end) + Bd_0*u_dt(i) + c_0];
    
end
% x_koop = Clift * xlift; % Koopman predictions
x_koop = X_koop;

%% ****************************  Plots  ***********************************
lw = 2;

figure;
subplot(2,1,1);
plot(x_true(1,:), x_true(2,:), 'k-', 'LineWidth', 1); hold on; grid on;
plot(x_koop(1,:), x_koop(2,:), 'r--', 'LineWidth', lw);
% plot(X_loc_x0(1,:), X_loc_x0(2,:), 'g--', 'linewidth', lw);
% plot(X_loc_0(1,:), X_loc_0(2,:), 'b--', 'linewidth', lw);
LEG = legend('True','Koopman','location','southeast');
set(LEG,'interpreter','latex')
title('Koopman Approximation of Dubin A/C North/East') 
xlim([min(x_true(1,:)),1.1*max(x_true(1,:))]);
ylim([min(x_true(2,:)), 1.1*max(x_true(2,:))]);

subplot(2,1,2);
lw = 1;
k_err = mod(x_koop(3,:) - mod(x_true(3,:), 2*pi), 2*pi);
plot(x_koop(3,:) - x_true(3,:), 'k-', 'LineWidth', lw); hold on; grid on;
% plot(X_loc_x0(3,:) - x_true(3,:), 'g--', 'LineWidth', lw); 
% plot(X_loc_0(3,:) - x_true(3,:), 'b--', 'LineWidth', lw); 
LEG = legend('Koopman','location','southeast');
set(LEG,'interpreter','latex')
title('Heading Angle Errors') 
% ylim([1.1*min(k_err), 1.1*max(k_err)]);

figure;
k_pos_error = vecnorm(x_true(1:2,:) - x_koop(1:2,:)); hold on; grid on;
kerr_sum = cumsum(k_pos_error);
plot(kerr_sum, 'k', 'LineWidth',2);
xlabel('Step')
ylabel('Cumulative Error');
title('Koopman Approximation Cumulative Position Error') 
