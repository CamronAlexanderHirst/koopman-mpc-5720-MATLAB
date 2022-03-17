clear all
close all
addpath('./Resources')
addpath('./Resources/qpOASES-3.1.0/interfaces/matlab') 
rng(2141444)

% Introduction:
% In this demo, we experiment with using the Koopman operator dynamics
% model of the aircraft for path-following control. We begin similarly to
% the DubinsAC.m demo, defining the dynamics, basis functions, and fitting
% an approximation of the Koopman operator. The Koopman operator is then
% used for model predictive control. Users can experiment with parameters
% regarding data collection, Koopman operator, and the MPC itself.
% Comparisons to a local linearization are presented.
% 
% !!! Note: due to the formulation of the MPC code, we must regress the
% position dependent dynamics. Future work would implement the position
% invariant model for MPC, but this would require a fairly significant
% re-write of the MPC code. Without the encoded model bias, we require more 
% data and basis functions to regress a decent model !!!

%% *************************** Exercises *******************************

% Data gathering params
Nsim = 500;  % number of data-gathering simulations
Ntraj = 500;  % number of steps in each data-gathering simulation

% Open-loop sim params:
T_max_pred = 5;  % Open-loop prediction simulation time (seconds)

% Closed-loop sim params:
Tpred = 3;  % MPC horizon (sec)
Tmax = 10;  % MPC closed-loop simlation time (seconds)
REF = 'line';  % Path-following reference ('line' or 'circle')
control_limits = 1;  % +- heading rate control constraint (rad/sec)

% Can also change the basis parameters like in DubinsAC.m!

% Questions:

% 0. Write down or save off all defaults. Initialize each question with the
% default parameters, then follow instructions.

% 1. Run the code. Look at the open-loop predictions in this example. Why 
% are the Koopman predictions worse than in DubinsAC.m? Why does MPC still
% kind of work?

% 2. Increase the MPC horizon (i.e. 5 sec). Decrease the MPC horizon i.e. 
% (1 sec) How does this change performance? Why?

% 3. Half Ntraj, then double it.How does changing the data gathering 
% parameters affect open-loop predictions? Closed-loop MPC performance?

% 4. How do solve times compare between Local linearized and Koopman MPC?
% Why the discrepancy?

% 5. Change the path-following reference to 'circle'. Why the discrepancy
% in the change in performance?

% 6. Change the control limits to a different value between 0-2. Why the
% change in performance?

% 7. Other things to try: mess with the basis functions. try commenting
% lines 99 thru 101 and uncommenting lines 103 thru 105. Look at the
% difference between these lifting functions. Why the change in
% performance?


%% *************************** Dynamics ***********************************

v = 1.;  % assumed constant velocity
% Assume zero control input
f_u =  @(t,x,u)([ v*cos(x(3,:)) ; v*sin(x(3,:)) ; u ] );
X0 = [0.; 0.; 0.];
n = 3; % number of states
m = 1; % number of control inputs

%% ************************** Discretization ******************************

deltaT = 0.1;  % step time, sec
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k3(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );


%% ************************** Basis functions *****************************

% Basis Function Definitions
n = 3;                                % Dimension of system (do not change)
% RBF 1 Centers
Nrbf = 100;                             % # of basis functions
cent = rand(n,Nrbf)*2 - 1;    % centers of each function
rbf_type = 'invquad'; 

% RBF 2 Centers
Nrbf2 = 0;                            % # of basis functions
cent2 = rand(n,Nrbf2)*2 - 1;          % centers of each function
rbf_type_2 = 'thinplate';             % type of function - one of 'thinplate', 'gauss', 'invquad', 'polyharmonic'

% Lifting mapping - RBFs + the state itself

extra_param = 1;
liftFun = @(xx)( [xx; rbf(xx,cent,rbf_type); rbf(xx, cent2,rbf_type_2, extra_param)]); 
Nlift = (Nrbf + Nrbf2) + n;

% extra_param = 1;
% liftFun = @(xx)( [xx; rbf(xx,cent,rbf_type); rbf(xx, cent2,rbf_type_2, extra_param); cos(xx(3,:)) ; sin(xx(3,:))]); 
% Nlift = (Nrbf + Nrbf2) + n + 2;


%% ************************** Collect data ********************************
tic
disp('Starting data collection')

% Random forcing
Ubig = 2*rand([Nsim Ntraj]) - 1;

% Random initial conditions
heading_sample_range = 2;     % +/- this range
Xcurrent = [rand(2,Ntraj)*2 - 1; rand(1,Ntraj)*heading_sample_range*2 - heading_sample_range];
X = []; Y = []; U = [];

for i = 1:Nsim
    Xnext = f_ud(0, Xcurrent, Ubig(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Xcurrent = Xnext;
end

fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ******************************* Lift ***********************************

disp('Starting LIFTING of dataset')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);
fprintf('Lifting dataset DONE, time = %1.2f s \n', toc);

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

Nsim = T_max_pred/deltaT;
u_dt = @(i)((-1).^(round(i/3))); % control signal

% Initial condition
x0 = [0.5; 0.5; 0.1];
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

% Local linearization predictor at 0
x = sym('x',[3;1]); u = sym('u',[1;1]);
Ac_0 = double(subs(jacobian(f_u(0,x,u),x),[x;u],[0;0;0;0]));
Bc_0 = double(subs(jacobian(f_u(0,x,u),u),[x;u],[0;0;0;0]));
c_0 = double(subs(f_u(0,x,u),[x;u],[0;0;0;0])) - Ac_0*[0;0;0] - Bc_0*0;
ABc = expm([Ac_0 [Bc_0 c_0] ; zeros(2,5)]*deltaT); % discretize
Ad_0 = ABc(1:3,1:3); Bd_0 = ABc(1:3,3); cd_0 = ABc(1:3,4); 
X_loc_0 = x0;


% Simulate
for i = 0:Nsim-1
    % Position-Variant Koopman Predictor
    xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)]; % Lifted dynamics
    
    % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
    
    % Local linearization predictor at x0
    X_loc_x0 = [X_loc_x0, Ad_x0*X_loc_x0(:,end) + Bd_x0*u_dt(i) + cd_x0];

    % Local linearization predictor at 0
    X_loc_0 = [X_loc_0, Ad_0*X_loc_0(:,end) + Bd_0*u_dt(i) + c_0]; % probably remove this
    
end
x_koop = Clift * xlift; % Koopman predictions Clift(1,:)



%% ****************************  Predictor Plots  *************************

lw = 2;
figure;
plot(x_true(1,:), x_true(2,:), 'k-', 'LineWidth', 1); hold on; grid on;
plot(x_koop(1,:), x_koop(2,:), 'r--', 'LineWidth', lw);
plot(X_loc_x0(1,:), X_loc_x0(2,:), 'g--', 'linewidth', lw);
plot(X_loc_0(1,:), X_loc_0(2,:), 'b--', 'linewidth', lw);
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','northeast');
set(LEG,'interpreter','latex')
xlim([-0.5,6]);
ylim([0.3, 1]);
xlabel('x');
ylabel('y');
title('Open-Loop Predictions of Dubin A/C Position');

%% ********************* Model Predictive Control *************************


Nsim = Tmax/deltaT;
switch REF
    case 'circle'
        radius = 10;
        dist = v*deltaT*Nsim / (2*pi*radius);
        
        yrr = radius * cos(2*pi*dist*[1:Nsim]/Nsim) - radius;  % y ref
        xrr = radius * sin(2*pi*dist*[1:Nsim]/Nsim);  % x ref
        ref_traj = [xrr; yrr; 0*[1:Nsim]];
        u_min = -1 * control_limits;  % control lower bound
        u_max = control_limits;  % control upper bound
        x0 = [0.0; 0.5; 0.0];  % initial state
    case 'line'
        yrr = 0*[1:Nsim];  % y ref
        xrr = v*deltaT*[1:Nsim] - 1;  % x ref
        ref_traj = [xrr; yrr; 0*[1:Nsim]];
        u_min = -1 * control_limits;  % control lower bound
        u_max = control_limits;  % control upper bound
        x0 = [-1.0; 2.0; 0.0];  % initial state
end


% Define Koopman controller
C = zeros(1,Nlift); C(1) = 1; C(2) = 1; C(3) = 0;  % define J wrt position
C = diag(C);
% Weight matrices
Q = diag([1,1,0]);  % state cost [x,y,heading]
R = 0.01;  % control cost
% Prediction horizon
Np = round(Tpred / deltaT);  % number of knots in MPC horizon

% Constraints
xlift_min = [nan(Nlift,1)];
xlift_max = [nan(Nlift,1)];

% Build Koopman MPC controller
koopmanMPC  = getMPC(Alift,Blift,Clift,0,Q,R,Q,Np, u_min, u_max, xlift_min, xlift_max,'qpoases');

x_koop = x0; x_loc = x0; x_koop_PI = x0;
XX_koop = x0; UU_koop = [];
XX_koop_PI = x0; UU_koop_PI = [];
XX_loc = x0; UU_loc = [];
time_koop = []; time_loc = [];

% Get Jacobian of the true dynamics (for local linearization MPC)
x = sym('x',[3 1]); syms u;
f_ud_sym = f_ud(0,x,u);
u_loc = 0;
Jx = jacobian(f_ud_sym,x);
Ju = jacobian(f_ud_sym,u);
Cy = diag([1 1 0]); % Output matrix: y = Cy*x
y_min = [nan(3,1)];
y_max = [nan(3,1)];


wasinfeas= 0;
ind_inf = [];

% Closed-loop simultion start
for i = 0:Nsim-1
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', i, Nsim)
    end
    
    % Current value of the reference signal
    yr = ref_traj(:,i+1);
    
    tic
    % Koopman MPC - Position variant
    xlift = liftFun(x_koop); % Lift
    u_koop = koopmanMPC(xlift,yr); % Get control input
    x_koop = f_ud(0,x_koop,u_koop); % Update true state
    time_koop = [time_koop toc];

    
    % Local linearization MPC
    tic
    Aloc = double(subs(Jx,[x;u],[x_loc;u_loc])); % Get local linearization
    Bloc = double(subs(Ju,[x;u],[x_loc;u_loc]));
    cloc = double(subs(f_ud_sym,[x;u],[x_loc;u_loc])) - Aloc*x_loc - Bloc*u_loc;
    [U_loc,~,optval] = solveMPCprob(Aloc,Bloc,Cy,cloc,Q,R,Q,Np,-1, 1,y_min,y_max,x_loc,yr); % Get control input
    u_loc = U_loc(1:m,1);
    if(optval == Inf) % Detect infeasibility
        ind_inf = [ind_inf i];
        wasinfeas = 1;
    end
    x_loc = f_ud(0,x_loc,u_loc); % Update true state
    time_loc = [time_loc toc];
        
    % Store values
    XX_koop = [XX_koop x_koop];
    UU_koop = [UU_koop u_koop];
    XX_loc = [XX_loc x_loc];
    UU_loc = [UU_loc u_loc];
end

if(isempty(ind_inf))
    ind_inf = Nsim;
end

%% ************************** MPC Plots **********************************
figure; 
plot(ref_traj(1,:), ref_traj(2,:), 'k--', 'LineWidth', lw); hold on; grid on;
plot(XX_koop(1,:), XX_koop(2,:), 'g-', 'LineWidth', lw);
plot(XX_loc(1,:), XX_loc(2,:), 'r--', 'LineWidth', lw);
LEG = legend('Ref Path','Koopman','Local at $x_0$','location','northeast');
set(LEG,'interpreter','latex');
xlabel('x');
ylabel('y');
title('Closed-Loop Dubin A/C Path-Following MPC');

figure; 
plot(time_koop(:), 'g-', 'LineWidth', lw);hold on; grid on;
plot(time_loc(:), 'r--', 'Linewidth', lw)
LEG = legend('Koopman','Local at $x_0$','location','northeast');
set(LEG,'interpreter','latex')
xlabel('simulation step');
ylabel('Solve Time (s)');
title('MPC Solve Times');

figure;
plot(UU_koop(:), 'g-', 'LineWidth', lw);hold on; grid on;
plot(UU_loc(:), 'r--', 'Linewidth', lw)
LEG = legend('Koopman','Local at $x_0$','location','northeast');
set(LEG,'interpreter','latex')
xlabel('simulation step');
ylabel('Yaw Rate Control Signal');
title('MPC Control Signal');