clear all
close all
addpath('./Resources')
rng(2141444)

%% TODO: Add some more commenting here!
%%%%%%%%%%%%%%%%%%%%%%%%%%% Data collection parameters

Nsim = 10;                    % # of simulations
Ntraj = 50;                   % # of steps in each simulation
heading_sample_range = 1;     % +/- this range

%%%%%%%%%%%%%%%%%%%%%%%%%%% Basis functions

% Problem -- try out other rbs, with 'invquad', 'gauss', 'polyharmonic' -- 
% Can also try combos of these -- can you find a basis that gives the
% best results?

% Basis Function Definitions
n = 3;                                % Dimension of system (do not change)
% RBF 1 Centers
Nrbf = 1;                             % # of basis functions
cent = rand(n,Nrbf)*2 - 1;            % centers of each function
rbf_type = 'invquad'; 

% RBF 2 Centers
Nrbf2 = 0;                            % # of basis functions
cent2 = rand(n,Nrbf2)*2 - 1;          % centers of each function
rbf_type_2 = 'thinplate';             % type of function - one of 'thinplate', 'gauss', 'invquad', 'polyharmonic'

% Lifting mapping - RBFs + the state itself
extra_param = 1;
liftFun = @(xx)( [xx; rbf(xx(3,:),cent,rbf_type); rbf(xx(3,:), cent2,rbf_type_2, extra_param);]);
Nlift = (Nrbf + Nrbf2) + n;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation Parameters
% Simulation Results
Tmax = 30;                   % Max Timestep

% Once you have a good approximation, try changing the intial heading to pi
% -- what happens? Can you fix it?
x0 = [5.; 2.; 0.];       % Initial Condition (North / East / Heading)  
        
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
k4 = @(t,x,u) ( f_u(t,x + k3(t,x,u)*deltaT,u) );        % k1 should be k3
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

% In MPC, compare NMPC w/ Koopman-based MPC
% Example with Dubin's

%% ************************** Collect data ********************************
tic
disp('Starting data collection')

% Random forcing
Ubig = 2*rand([Nsim Ntraj]) - 1;

% Random initial heading
Xcurrent = [zeros(2,Ntraj); rand(1,Ntraj)*heading_sample_range*2 - heading_sample_range];


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
u_dt = @(i)((-1).^(round(i/250))); % control signal

% u_dt = @(i)(1.);

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
% Simulate
for i = 0:Nsim-1
    x_tmp = [0; 0; X_koop(3,end)];
    x_tmp_lift = liftFun(x_tmp);
    delta_x_lift = Alift*x_tmp_lift + Blift*u_dt(i);
    X_koop = [X_koop [X_koop(1:2,end)+delta_x_lift(1:2); delta_x_lift(3)]];
    
 % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
    
%     % Local linearization predictor at x0
    X_loc_x0 = [X_loc_x0, Ad_x0*X_loc_x0(:,end) + Bd_x0*u_dt(i) + cd_x0];
%     
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
plot(X_loc_x0(1,:), X_loc_x0(2,:), 'g--', 'linewidth', lw);
plot(X_loc_0(1,:), X_loc_0(2,:), 'b--', 'linewidth', lw);
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southeast');
set(LEG,'interpreter','latex')
title('Koopman Approximation of Dubin A/C North/East') 
xlim([0.9*min(x_true(1,:)),1.1*max(x_true(1,:))]);
ylim([0.9*min(x_true(2,:)), 1.1*max(x_true(2,:))]);

subplot(2,1,2);
lw = 1;
k_err = x_koop(3,:) - x_true(3,:);
plot(x_koop(3,:) - x_true(3,:), 'k-', 'LineWidth', lw); hold on; grid on;
plot(X_loc_x0(3,:) - x_true(3,:), 'g--', 'LineWidth', lw); 
plot(X_loc_0(3,:) - x_true(3,:), 'b--', 'LineWidth', lw); 
LEG = legend('Koopman','Local at $x_0$','Local at 0','location','southeast');
set(LEG,'interpreter','latex')
title('Heading Angle Errors') 
ylim([1.1*min(k_err), 1.1*max(k_err)]);

figure;
k_pos_error = vecnorm(x_true(1:2,:) - x_koop(1:2,:)); hold on; grid on;
kerr_sum = cumsum(k_pos_error);
plot(kerr_sum, 'k', 'LineWidth',2);
xlabel('Step')
ylabel('Cumulative Error');
title('Koopman Approximation Cumulative Position Error') 
