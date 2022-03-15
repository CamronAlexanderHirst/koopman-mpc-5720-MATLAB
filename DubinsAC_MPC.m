clear all
close all
addpath('./Resources')
addpath('./Resources/qpOASES-3.1.0/interfaces/matlab') 
rng(2141444)


%% *************************** Dynamics ***********************************

v = 1.;  % assumed constant velocity

% Assume zero control input
f_u =  @(t,x,u)([ v*cos(x(3,:)) ; v*sin(x(3,:)) ; u ] );

X0 = [0.; 0.; 0.];
n = 3; % number of states
m = 1; % number of control inputs

%% ************************** Discretization ******************************

deltaT = 0.01;  % step time, sec

%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k3(t,x,u)*deltaT,u) );        % k1 should be k3
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

% In MPC, compare NMPC w/ Koopman-based MPC
%% ************************** Basis functions *****************************

basisFunction = 'rbf';
% RBF centers
Nrbf = 200;                            % Radial basis function -- # of bases
cent = rand(n,Nrbf)*2 - 1;
rbf_type = 'thinplate'; 
% Lifting mapping - RBFs + the state itself
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
Nlift = Nrbf + n;


%% ************************** Collect data ********************************
tic
disp('Starting data collection')
Nsim = 200;
Ntraj = 1000;

% Random forcing
Ubig = 2*rand([Nsim Ntraj]) - 1;

% Random initial conditions
Xcurrent = (rand(n,Ntraj)*2 - 1);
X = []; Y = []; U = [];

for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Ubig(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Xcurrent = Xnext;
end

% n_zeta = (nD+1)*ny + nD*m; % dimension of the delay-embedded "state"
% for i = 1:Nsim
%     Xnext = f_ud(0,Xcurrent,Ubig(i,:));
%     zeta_prev = zeta_current;
%     zeta_current = [[Cy*Xnext ; Ubig(i,:)] ; zeta_current( 1:end-ny-m , : ) ];
%     if(i > nD)
%         X = [X zeta_prev];
%         Y = [Y zeta_current];
%         U = [U Ubig(i,:)];
%     end
%     Xcurrent = Xnext;
% end
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

Tmax = 7;
Nsim = Tmax/deltaT;
u_dt = @(i)((-1).^(round(i/30))); % control signal

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

% Delayed initial condition (assume random control input in the past)
% x = x0;
% xstart = [Cy*x ; NaN(nD*(ny+m),1)];
% for i = 1:nD
%     urand = 2*rand(m,1) - 1;
%     xp = f_ud(0,x,urand);
%     xstart = [Cy*xp ; urand; xstart(1:end-ny-m)];
%     x = xp;
% end
% xlift = liftFun(xstart);


% Simulate
for i = 0:Nsim-1
    % Koopman predictor
    xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)]; % Lifted dynamics
    
    % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
    
    % Local linearization predictor at x0
    X_loc_x0 = [X_loc_x0, Ad_x0*X_loc_x0(:,end) + Bd_x0*u_dt(i) + cd_x0];     
    % Local linearization predictor at 0
    X_loc_0 = [X_loc_0, Ad_0*X_loc_0(:,end) + Bd_0*u_dt(i) + c_0];
    
end
x_koop = Clift * xlift; % Koopman predictions Clift(1,:)



%% ****************************  Predictor Plots  *************************

lw = 2;
figure
title('Koopman Approximation of Dubin A/C')
plot(x_true(1,:), x_true(2,:), 'k-', 'LineWidth', 1); hold on; grid on;
plot(x_koop(1,:), x_koop(2,:), 'r--', 'LineWidth', lw);
plot(X_loc_x0(1,:), X_loc_x0(2,:), 'g--', 'linewidth', lw);
plot(X_loc_0(1,:), X_loc_0(2,:), 'b--', 'linewidth', lw);
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','northeast');
set(LEG,'interpreter','latex')
xlim([-0.5,6]);
ylim([0.3, 1]);

%% ********************* Model Predictive Control *************************
% disp('Press any key for model predictive control')
% pause

Tmax = 5; % Simlation length (seconds)
Nsim = Tmax/deltaT;
REF = 'circle'; 
switch REF
    case 'circle'
        radius = 3;
        dist = v*deltaT*Nsim / (2*pi*radius);
        
        yrr = radius * cos(2*pi*dist*[1:Nsim]/Nsim);  % y ref
        xrr = radius * sin(2*pi*dist*[1:Nsim]/Nsim);  % x ref
        ref_traj = [xrr; yrr; 0*[1:Nsim]];
        u_min = -0.5;  % control lower bound
        u_max = 0.5;  % control upper bound
        x0 = [0.0; 3.5; 0.0];  % initial state
    case 'line'
        radius = 3;
        yrr = 0*[1:Nsim];  % y ref
        xrr = v*deltaT*[1:Nsim];  % x ref
        ref_traj = [xrr; yrr; 0*[1:Nsim]];
        u_min = -0.5;  % control lower bound
        u_max = 0.5;  % control upper bound
        x0 = [0.0; 2.0; 0.2];  % initial state
end


% Define Koopman controller
C = zeros(1,Nlift); C(1) = 1; C(2) = 1; C(3) = 0;  % define J wrt position
C = diag(C);
% Weight matrices
Q = diag([1,1,0]);  % state cost
%Q = 1;
R = 0.01;  % control cost
% Prediction horizon
Tpred = 1;  % MPC horizon (sec)
Np = round(Tpred / deltaT);  % number of knots in MPC horizon

% Constraints
xlift_min = [nan(Nlift,1)];
xlift_max = [nan(Nlift,1)];

% Build Koopman MPC controller
koopmanMPC  = getMPC(Alift,Blift,Clift,0,Q,R,Q,Np, u_min, u_max, xlift_min, xlift_max,'qpoases');

x_koop = x0; x_loc = x0;
XX_koop = x0; UU_koop = [];
XX_loc = x0; UU_loc = [];

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
    
    % Koopman MPC
    xlift = liftFun(x_koop); % Lift
    u_koop = koopmanMPC(xlift,yr); % Get control input
    x_koop = f_ud(0,x_koop,u_koop); % Update true state
    
    % Local linearization MPC
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
    
    % Nonlinear MPC
    
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
figure; title('Koopman Approximation of Dubin A/C')
plot(ref_traj(1,:), ref_traj(2,:), 'k--', 'LineWidth', lw); hold on; grid on;
plot(XX_koop(1,:), XX_koop(2,:), 'g-', 'LineWidth', lw); 
plot(XX_loc(1,:), XX_loc(2,:), 'r--', 'LineWidth', lw);
LEG = legend('ref','Koopman','Local at $x_0$','location','northeast');
set(LEG,'interpreter','latex')
% xlim([-0.5,6]);
% ylim([0.3, 1]);