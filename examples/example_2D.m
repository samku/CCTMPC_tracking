addpath(genpath('../src'))

% ---------------
% Note: we require that those softwares are already installed on the system
addpath(genpath('~/workSW/gurobi1100/linux64/matlab/'))
run('~/Documenti/CCTMPC/tbxmanager/startup.m')
addpath('~/workSW/casadi-3.6.4-linux64-matlab2018b/')
import casadi.*
% ---------------


% Dynamics
A=(1/5)*[5 1; -1 4];    nx = size(A,2);
B=(1/5)*[0;1];          nu = size(B,2);
Bw=(1/5)*[1 0; 0 1];    nw = size(Bw,2);

% Constraints
HX = [eye(nx); -eye(nx)]; hX = [6.8;10; 10;5];  mX = size(HX,1);
HU = [eye(nu); -eye(nu)]; hU = [10; 10];        mU = size(HU,1);
HW = [eye(nx); -eye(nx)]; hW = [0.5;2; 0.5;2];  mW = size(HW,1);

% Output
C = [1, 0]; ny = size(C,1);
D = zeros(ny, nu);

% LTI system
A_convh = {}; B_convh = {};
A_convh{end+1} = A; B_convh{end+1} = B;

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);

% we carry this representation for disturbance, as the number of vertices
% is <= compared to the original one (it's a transformation matrix)
W_dist = Bw*Polyhedron(HW, hW); 

dynSystem = DynSystem(A_convh, B_convh, C, D, X, U, W_dist);

% get Configuration-Constrained polytope
C_tilde = getConfConstrTemplate(nx,4);
ccpoly = CCPolytope(dynSystem, C_tilde, C_tilde);

% Cost Function definition: we can pass either structs of matrices, or
% function handle directly defined (just be sure that the function is
% -strictly- convex, and that the function signatures are:
% OTI_cost : @(ys, us, r, th_s)
% Stage/Term_cost :  @(y_k,u_k,ys,us)
% TODO!!!: decide if us is a matrix or a vector!
RCI_cost = struct(  "Qv",blkdiag(10*eye(nx),eye(nu)),...
    "Qc",blkdiag(eye(nx),eye(nu)),...
    "Qr",eye(ny) );
% Qcost = eye(nx);
% gamma = 0.95;
% Pcost = Qcost*((1-gamma^2)^-1) ; % or better separate the two costs? or use a varargin?

gamma = 0.95;
[Qcost_han, Pcost_han] = defineCustomTrackCost(dynSystem, ccpoly, gamma);

% TODO: should we consider u_k as the paper? or as a 2D matrix, which expands on the 3rd dimension?
costFunMan = CostFunctionManager(dynSystem, ccpoly, RCI_cost, Qcost_han, Pcost_han); % is it needed? in theory to be general yes

% build the MPC scheme (under the hood, it instantiate an OptimalController)

solver_opts = {'solver','gurobi','verbose',0};
optController = OptimalController(dynSystem, ccpoly, costFunMan, 10, gamma, solver_opts);
cctmpc = CCTMPC;

N_mpc = 100;

% get Initial state x0

% TODO: build a reference generator

for i=1:N_mpc
    % simulate the system. Also, call random disturbances (to be added as a
    % method for the system).
    
    % take care of specific needs, i.e.: changing references at each time
    % step; update model at each time step, and so on
    
end

function C_tilde = getConfConstrTemplate(nx,m)
% this ensure to get minimal representation, m_bar = m
% by using polar coordinates
assert(nx==2)
C_tilde = zeros(m,nx);
for i=1:m
    phi = 2*(i-1)*pi/m;
    C_tilde(i,:) = [cos(phi), sin(phi)];
end
end

function [Qtrack_han, Ptrack_han] = defineCustomTrackCost(sys,ccpoly,gamma)
weightSq2Norm = @(vector, mat) vector'*mat*vector;

Qx_dist=10*eye(sys.nx);
Qu_dist=eye(sys.nu);

V_bar = mean(cat(3, ccpoly.Vi_s{:}), 3); % each Vi_s{i} is a 2D matrix, concatenating on 3rd dim.
U_bar = (1/ccpoly.m_bar)*repmat(eye(sys.nu),1,ccpoly.m_bar); % nu x (nu*m_bar)

Qy_dist=zeros(ccpoly.m,ccpoly.m);
for i=1:m_bar
    Qy_dist=Qy_dist+(ccpoly.Vi_s{i}-V_bar)'*Qx_dist*(ccpoly.Vi_s{i}-V_bar);
end

Qub_dist=zeros(sys.nu*ccpoly.m_bar,sys.nu*ccpoly.m_bar);
Iu_mats={};
for i=1:ccpoly.m_bar
    Iu_mats{i}=sparse(sys.nu,sys.nu*ccpoly.m_bar);
    Iu_mats{i}(:,(i-1)*sys.nu+1:i*sys.nu)=eye(sys.nu);
    Qub_dist=Qub_dist+(Iu_mats{i}-U_bar)'*Qu_dist*(Iu_mats{i}-U_bar);
end

Q=blkdiag(Qy_dist,Qub_dist)+0.005*eye(m+m_bar*nu,m+m_bar*nu);
P = (1/(1-gamma^2))*Q;

Qtrack_han = @(y_k,u_k,ys,us) weightSq2Norm([y_k-ys; u_k(:)-us(:)],Q);
Ptrack_han = @(y_N,u_N,ys,us) weightSq2Norm([y_N-ys; u_N(:)-us(:)],P);

end




