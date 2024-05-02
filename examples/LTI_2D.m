addpath(genpath('../src'))

% -----------------------
% NOTE: we require that those softwares are already installed on the system
addpath(genpath('~/workSW/gurobi1100/linux64/matlab/'))
% NOTE : at least cplex and mosek are required
addpath(genpath('~/workSW/mosek/10.1/toolbox/r2017aom/'))
addpath(genpath('/opt/ibm/ILOG/CPLEX_Studio1210/cplex/matlab/x86-64_linux'))
run('~/Documenti/CCTMPC/tbxmanager/startup.m')
addpath('~/workSW/casadi-3.6.4-linux64-matlab2018b/')
import casadi.*
% -----------------------

%% System definition
% Dynamics
A = [1.1 1; 0 1];   nx = size(A,2);
B = [0.5;1];        nu = size(B,2);
Bw = [1 0; 0 1];    nw = size(Bw,2);

% Constraints
HX = [eye(nx); -eye(nx)]; hX = [5;3;5;2];       mX = size(HX,1);
HU = [eye(nu); -eye(nu)]; hU = [2;1];           mU = size(HU,1);
HW = [eye(nx); -eye(nx)]; hW = [0;0.5;0;0.5];   mW = size(HW,1);

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

% since it's an LTI system, we ignore additional parameters
A_curr_han = @(A_convh,params) A_convh{1};
B_curr_han = @(B_convh,params) B_convh{1};

sys = DynSystem(A_convh, B_convh, C, D, X, U, W_dist,A_curr_han,B_curr_han);
% -----------------------

% By using polar coordinates for parameterization, we ensure to get
% a CCPolytope with minimal representation, i.e. m_bar = m
C_tilde = getTemplateMatrix(nx,16);
ccPoly = CCPolytope(sys, C_tilde);

% Cost Function definition: we can pass either structs of matrices, or
% function handles defined externally. Just be sure that the function is
% -strictly- convex, and that the function signatures are:
% RCI_cost : @(ys, us, r, th_s)
% Stage/Term_cost : @(y_k/N,u_k/N,ys,us)
RCI_cost = struct(  "Qv",blkdiag(10*eye(nx),eye(nu)),...
    "Qc",blkdiag(eye(nx),eye(nu)),...
    "Qr",eye(ny) );

gamma = 0.95;
[Qcost_han, Pcost_han] = defineCustomTrackCost(sys, ccPoly, gamma);

% lump all cost function definitions together
costFunMan = CostFunctionManager(sys, ccPoly, RCI_cost, Qcost_han, Pcost_han);

% build the MPC scheme
solver_opts = {'solver','gurobi','verbose',0};
N_ocp = 5;
variableConvh = false; % used as flag for OptimalController class (faster code)
useBasicOCP = false;
cctmpc = CCTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, solver_opts, variableConvh,useBasicOCP);

%% Define the control schemes for comparison

lambda = 0.999; % used for compute a tight approx. of the mRPI set for terminal set.

elastic_mpc = ETMPC_tracking(sys, ccPoly, RCI_cost, N_ocp, lambda, solver_opts);
limon_mpc = RTMPC_Limon(sys, ccPoly, RCI_cost, N_ocp, lambda, solver_opts);

%% Simulation

N_mpc = 25+1; % simulation steps

% (state space) system dynamics
x_sys = zeros(nx, N_mpc+1); x_sys(:,1) = [-3.43;3]; % initial condition
u_sys = zeros(nu, N_mpc);
w_sys = zeros(nx, N_mpc);

% define reference to be followed
r_sys = zeros(ny, N_mpc); r_change = 12+1;
r_sys(:,1:r_change-1) = 5; r_sys(:,r_change:end) = -5;

% (parameter space) optimal control dynamics
OCP_y = cell(1,N_mpc); OCP_u = cell(1,N_mpc);
OCP_ys = cell(1,N_mpc); OCP_us = cell(1,N_mpc);

RCI_ys = cell(1,N_mpc); RCI_us = cell(1,N_mpc);

Lyap_cost = zeros(1,N_mpc);


for t = 1:N_mpc
    disp("MPC iteration: " + t+"/"+N_mpc)
    
    % compute input for the system
    u_sys(:,t) = cctmpc.solve(x_sys(:,t), r_sys(:,t));
    
    % compute worst case disturbance (with respect to the origin)
    w_sys(:,t) = getWorstCaseDist(sys, x_sys(:,t), u_sys(:,t));
    
    % evolve the system
    x_sys(:,t+1) = sys.A_curr()*x_sys(:,t) + sys.B_curr()*u_sys(:,t) + w_sys(:,t);
    
    % save computed data
    [OCP_y{t}, OCP_u{t}, OCP_ys{t}, OCP_us{t}] = cctmpc.ocpSol{:};
    [RCI_ys{t}, RCI_us{t}] = cctmpc.rciSol{:};
    Lyap_cost(t) = cctmpc.Lyapunov_cost;
end

%% compute the feasibility regions, and compare it with the Maximal Robust
% Control Invariant set computed according to Rakovic approach.
feasRegionCCTMPC = cctmpc.computeFeasRegion(200);
feasRegionETMPC = elastic_mpc.computeFeasRegion(200);
feasRegionRTMPC = limon_mpc.computeFeasRegion(200);

MRCI = computeMRCIset(sys);

% compare feasible region distance from the MRCI
dist_CCTMPC_MRCI = setDistance(feasRegionCCTMPC, MRCI);
dist_ETMPC_MRCI = setDistance(feasRegionETMPC, MRCI);
dist_RTMPC_MRCI = setDistance(feasRegionRTMPC, MRCI);

disp("Hausdorff distance (MRCI,CCTMPC): " + dist_CCTMPC_MRCI)
disp("Hausdorff distance (MRCI,ETMPC_tracking): " + dist_ETMPC_MRCI)
disp("Hausdorff distance (MRCI,Limon_RTMPC): " + dist_RTMPC_MRCI)

%% Plot section
% defining the plot scripts separately to improve code readibility
LTI_2D_PlotPhasePlane;
LTI_2D_PlotLyapunov;

%% ------------------------- support functions ----------------------------

function [Qtrack_han, Ptrack_han] = defineCustomTrackCost(sys,ccPoly,gamma)
% We provide an instance of a user-specified cost function, passed as an
% handle. In this case, we seek a tradeoff between tracking performance and
% CCPolytope volume. See "2D Example" in the Section IV from the paper.

weightSq2Norm = @(vector, mat) vector'*mat*vector;
nx = sys.nx; nu = sys.nu;
m = ccPoly.m; m_bar = ccPoly.m_bar;

% RCI cost matrices (in state space)
Qv_x = 10*eye(nx);
Qv_u = eye(nu);

V_bar = mean(cat(3, ccPoly.Vi_s{:}), 3); % each Vi_s{i} is a 2D matrix, averaging along 3rd dim.
U_bar = (1/m_bar)*repmat(eye(nu),1,m_bar); % nu x (nu*m_bar)

% cumulating matrix cost (in parameter space)
Qv_yp = zeros(m,m);
for i = 1:m_bar
    Qv_yp = Qv_yp+(ccPoly.Vi_s{i}-V_bar)'*Qv_x*(ccPoly.Vi_s{i}-V_bar);
end

Qv_up = zeros(nu*m_bar, nu*m_bar);
Iu_mats = {};
for i = 1:m_bar
    Iu_mats{i} = sparse(nu,nu*m_bar);
    Iu_mats{i}(:,(i-1)*nu+1:i*nu) = eye(nu);
    Qv_up = Qv_up+(Iu_mats{i}-U_bar)'*Qv_u*(Iu_mats{i}-U_bar);
end

Q = blkdiag(Qv_yp,Qv_up)+0.001*eye(m+m_bar*nu,m+m_bar*nu);
P = (1/(1-gamma^2))*Q;

% define function handles (here quadratic costs)
Qtrack_han = @(y_k,u_k,ys,us) weightSq2Norm([y_k-ys; u_k(:)-us(:)],Q);
Ptrack_han = @(y_N,u_N,ys,us) weightSq2Norm([y_N-ys; u_N(:)-us(:)],P);

end


function w_dist = getWorstCaseDist(sys,x_curr,uMPC)
W_dist_vert = sys.W_dist.V';

x_next = [];
dist_to_sp = [];
for j=1:size(W_dist_vert,2) % enumerate disturbance realizations
    x_next(:,j) = sys.A_curr*x_curr + sys.B_curr*uMPC + W_dist_vert(:,j);
    dist_to_sp(j) = norm(x_next(:,j),2)^2;
end
[~,id_worst] = max(dist_to_sp);
w_dist = W_dist_vert(:,id_worst);
end


function MRCI = computeMRCIset(sys)
% compute the Maximal Robust Control Invariant set, by starting from the
% state constraint set and taking into account system dynamics.
MRCI = sys.X;

for i=1:100 % generally a lot less iterations are required
    d = sys.W_dist.support(MRCI.A');
    
    H_LHS = [MRCI.A*sys.A_curr(), MRCI.A*sys.B_curr();
        sys.X.A, zeros(size(sys.X.A,1),sys.nu);
        zeros(size(sys.U.A,1),sys.nx), sys.U.A];
    h_RHS = [MRCI.b-d;
        sys.X.b;
        sys.U.b];
    MRCI_next = Polyhedron(H_LHS,h_RHS).projection(1:sys.nx).minHRep().minVRep();
    
    if MRCI_next == MRCI % we converged to the actual set.
        break
    else
        MRCI = MRCI_next;
    end
end

end


function dist = setDistance(polytope1, polytope2)
% compute the l-norm Hausdorff distance between polytopes
vert_poly2 = polytope2.V'; % (nx,m_bar)
[nx, m_bar] = size(vert_poly2);

x_pts = sdpvar(nx, m_bar,'full');
unit_ball_pts = sdpvar(nx, m_bar,'full');
eps = sdpvar(1,1,'full');

constr = [vert_poly2 == x_pts + unit_ball_pts;
    polytope1.A*x_pts <= repmat(polytope1.b, 1, m_bar);
    [eye(nx); -eye(nx)]*unit_ball_pts <= repmat(eps*ones(2*nx,1),1,m_bar)];

optimize(constr, eps, sdpsettings('verbose',0));
dist = value(eps);

end