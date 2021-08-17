%------------------------------------------%
% Object MPC
%------------------------------------------%

%% Dimensions

nx=29;  % No. of differential states
nu=16;  % No. of controls
nz=0;  % No. of algebraic states
ny=45; % No. of outputs
nyN=29; % No. of outputs at the terminal point
np=78; % No. of model parameters
nc=0; % No. of general constraints
ncN=0; % No. of general constraints at the terminal point
nbx = 29; % No. of bounds on states
nbu = 16; % No. of bounds on controls

% state and control bounds
nbx_idx = 1:29; % indexs of states which are bounded
nbu_idx = 1:16; % indexs of controls which are bounded

%% create variables
addpath('/home/akash/Downloads/casadi-linux-matlabR2014b-v3.5.5')
addpath('./kinematics/')
addpath('./kinematics/screw/screws/') %required for ad()
addpath('./kinematics/screw/util/') %required for isequalf()

load_init_params

o.reset();
Fc_hat = blkdiag(o.contacts(1).x_hat_, o.contacts(2).x_hat_, o.contacts(3).x_hat_, o.contacts(4).x_hat_);

% isunix = 1;
import casadi.*

states   = SX.sym('states',nx,1);   % differential states
controls = SX.sym('controls',nu,1); % control input
alg      = SX.sym('alg',nz,1);      % algebraic states
params   = SX.sym('paras',np,1);    % parameters
refs     = SX.sym('refs',ny,1);     % references of the first N stages
refN     = SX.sym('refs',nyN,1);    % reference of the last stage
Q        = SX.sym('Q',ny,1);        % weighting matrix of the first N stages
QN       = SX.sym('QN',nyN,1);      % weighting matrix of the last stage
aux      = SX.sym('aux',ny,1);      % auxilary variable
auxN     = SX.sym('auxN',nyN,1);    % auxilary variable

%% Dynamics

pos=states(1:3);
orien=states(4:6);
vel=states(7:9);
omega=states(10:12);  
gr=states(13);
lam=states(14:29);
u=controls;
Msym = reshape(params(1:36),6,6);
Csym = reshape(params(37:72),6,6);
Nsym = reshape(params(73:78),6,1);

% explicit ODE RHS
x_dot = [zyx2R(states(4:6))*states(7:9); zyx2E(states(4:6))*states(10:12);...
       Msym\(-Csym*states(7:12)-Nsym+o.G*Fc_hat*states(14:29)); 0; u];

% algebraic function
z_fun = [];                   

% implicit ODE: impl_f = 0
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;        
     
%% Objectives and constraints

% inner objectives
h = [pos;orien;vel;omega;gr;lam;u];
hN = h(1:nyN);

% outer objectives
obji = 0.5*(h-refs)'*diag(Q)*(h-refs);
objN = 0.5*(hN-refN)'*diag(QN)*(hN-refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);

% general inequality constraints
general_con = [];
general_con_N = [];

%% NMPC discretizing time length [s]

Ts_st = 0.01; % shooting interval time
