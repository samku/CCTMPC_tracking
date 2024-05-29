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
% Define the discrete LPV system
Ts = 0.025;
[A_is, B_is, Bw] = defDTLateralDynamics(Ts);

vx_interval  = [40,65]/3.6; % longitudinal speed [m/s]
[A_convh, B_convh, P_hat] = defPolytopicLPVEnclosure(A_is, B_is, vx_interval);

nx = size(A_convh{1},2);
nu = size(B_convh{1},2);
nw = size(Bw,2);

% Output
C = [1,0,0,0]; ny = size(C,1);
D = zeros(ny, nu);

% Constraints
X = Polyhedron([eye(nx); -eye(nx)], repmat([3; 4; 20*pi/180; 3],2,1) );
U = Polyhedron([eye(nu); -eye(nu)], repmat([10*pi/180; 2],2,1) );

W_dist = Ts*Bw*Polyhedron([eye(nw); -eye(nw)], [(10)^2;0]);

% LPV system representation (enclosed in the convex hull above defined)
A_curr_han = @(A_convh,vx) vx*A_is{1} + (1/vx)*A_is{2} + A_is{3};
B_curr_han = @(B_convh,vx) vx*B_is{1} + (1/vx)*B_is{2} + B_is{3};

sys = DynSystem(A_convh, B_convh, C, D, X, U, W_dist,A_curr_han,B_curr_han);

% compute a Configuration-Constrained polytope
C_tilde = [eye(nx); -eye(nx)];
ccPoly = CCPolytope(sys, C_tilde, C_tilde); % (sys,C_tilde,SetDistance)

% Cost Function definition: we can either pass structs of matrices, or
% function handles defined externally. We only require that the cost
% functions are -strictly- convex, and that the function signatures are:
% RCI_cost : @(ys, us, r, th_s)
% Stage/Term_cost : @(y_k/N,u_k/N,ys,us)
RCI_cost = struct(  "Qv",blkdiag(0.001*eye(nx), 0.01*eye(nu)),...
    "Qc",blkdiag(100*eye(nx), 0.01*eye(nu)),...
    "Qr",10*eye(ny) );

gamma = 0.95;
Qmat = blkdiag(0.001*eye(ccPoly.m),  0.01*eye(nu*ccPoly.v));
Pmat = ((1-gamma^2)^-1)*Qmat;

% collect all cost function definitions in a single class
costFunMan = CostFunctionManager(sys, ccPoly, RCI_cost, Qmat, Pmat);

% build the MPC scheme
solver_opts = {'solver','gurobi','verbose',0};
N_ocp = 3;
variableConvh = true; % flag for OptimalController class (fixed matrices produces faster code)
useBasicOCP = false;
cctmpc = CCTMPC(sys, ccPoly, costFunMan, N_ocp, gamma, solver_opts, variableConvh,useBasicOCP);

%% Simulation

N_mpc = 120+1; % simulation steps

% (state space) system dynamics
x_sys = zeros(nx, N_mpc+1); x_sys(:,1) = [-1.5;0;0;0]; % initial condition
u_sys = zeros(nu, N_mpc);
w_sys = zeros(nx, N_mpc);

% define reference to be followed
r_sys = zeros(ny, N_mpc); r_change = 30+1;
r_sys(:,1:r_change-1) = -1.5; r_sys(:,r_change:end) = 1.5;

% (parameter space) optimal control dynamics
OCP_y = cell(1,N_mpc); OCP_u = cell(1,N_mpc);
OCP_ys = cell(1,N_mpc); OCP_us = cell(1,N_mpc);

RCI_ys = cell(1,N_mpc); RCI_us = cell(1,N_mpc);

Lyap_cost = zeros(1,N_mpc);

% define longitudinal dynamic profile
longDyn = defLongDyn(Ts, vx_interval, N_mpc, r_change);

% main loop
for t = 1:N_mpc
    disp("MPC iteration: " + t+"/"+N_mpc)    
    % compute the current A/B_convh
    [A_convh_curr, B_convh_curr] = updateDynModel(A_is, B_is, P_hat, longDyn.vx(t));
    
    % compute input for the system
    u_sys(:,t) = cctmpc.solve(x_sys(:,t), r_sys(:,t),A_convh_curr, B_convh_curr);
    
    % compute random disturbance
    w_sys(:,t) = W_dist.randomPoint();
    
    % propagate system dynamics
    x_sys(:,t+1) = sys.A_curr(longDyn.vx(t))*x_sys(:,t) + sys.B_curr(longDyn.vx(t))*u_sys(:,t) + w_sys(:,t);
    
    % save computed data
    [OCP_y{t}, OCP_u{t}, OCP_ys{t}, OCP_us{t}] = cctmpc.ocpSol{:};
    [RCI_ys{t}, RCI_us{t}] = cctmpc.rciSol{:};
    Lyap_cost(t) = cctmpc.Lyapunov_cost;
end

%% Define polytope projections
projVerts = cell(1,6);
for t = 1:N_mpc
    % (y_cells,projs), each cell has two vertices
    ccPolyProjs = projPolytope(ccPoly.F,{OCP_y{t},OCP_ys{t},RCI_ys{t}},[1,3]);
    % iterate over time
    for i = 1:numel(ccPolyProjs)
        projVerts{i} = [projVerts{i}, ccPolyProjs{i}'];
    end
end
projVerts = cell2struct(projVerts,...
    {'ey_y','ey_ys','ey_rci','epsi_y','epsi_ys','epsi_rci'},2);

%% Plot section (Driving scenario)
% generate plots in separate scripts to improve code readibility
LateralDynamics_PlotStateSpace;
LateralDynamics_PlotLyapunov;

%% ------------------------- Support functions ----------------------------

function [A_is, B_is, Bw] = defDTLateralDynamics(Ts)
% Continuous time system definition with empirical parameters
param1=-171.2893;
param2=85.2523;
param3=42.1875;
param4=-199.6488;

param5=65.8919;
param6=43.6411;
param7=0.2287;

param8=0.0018;
param9=-0.0022;

A1=[0 0 1 0; 0 0 0 -1; 0 0 0 0; 0 0 0 0];
A2=[0 0 0 0; 0 param1 0 param2; 0 0 0 0; 0 param3 0 param4];
A0=[0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
nx=size(A0,2);

B=[0 0; param5 0; 0 0; param6 param7*10^(-3)];
nu=size(B,2);
Bw=[0;param8;0;param9];

% Euler discretization
A1_d = Ts*A1; A2_d = Ts*A2; A3_d = Ts*A0 + eye(nx);
B1_d = zeros(nx,nu); B2_d = B1_d; B3_d = Ts*B;

A_is = {A1_d, A2_d, A3_d};
B_is = {B1_d, B2_d, B3_d};
end


function [A_convh, B_convh, P_hat] = defPolytopicLPVEnclosure(A_is, B_is, vx_interval)
% compute polytopic over-approximation of the parameter set
p1_val = vx_interval;
p2_val = 1./p1_val; % 1/vx

p_vals = [p1_val; p2_val];

slope = (p_vals(2,1)-p_vals(2,2)) / (p_vals(1,1)-p_vals(1,2));

% encapsulate the function p2 = f(p1)
h = sdpvar(1); p1 = sdpvar(1); p2 = sdpvar(1);

cost = h;
constr = [p1*p2 == 1;
    p2 == slope*p1+h;
    p1_val(1) <= p1 <= p1_val(end)];
optimize(constr, cost, sdpsettings('verbose',0))

h = value(h);

% get the 4 vertices
v_lh = p_vals(:,1); v_hl = p_vals(:,2);
v_ll = [p_vals(1,1); slope*p_vals(1,1)+h];
v_hh = [(slope)\(p_vals(2,2)-h); p_vals(2,2)];

P_hat = Polyhedron([v_lh'; v_hl'; v_ll'; v_hh']);
P_hat_V = num2cell(P_hat.V',1); np = length(P_hat_V);

A_convh = cell(np,1); B_convh = cell(np,1);
for v=1:np
    p = [P_hat_V{v};1];
    Atmp = p(1)*A_is{1}+p(2)*A_is{2}+p(3)*A_is{3};
    Btmp = p(1)*B_is{1}+p(2)*B_is{2}+p(3)*B_is{3};
    A_convh{v} = Atmp;
    B_convh{v} = Btmp;
end

end

function vert_projs = projPolytope(F,y_cells,projs)
% project the CCPolytope for relevant state components
nx = size(F,2);
vert_projs = cell(length(y_cells),length(projs));

for p = 1:length(projs)
    for y=1:length(y_cells)
        poly_t = Polyhedron(F,y_cells{y}(:,1)); % use y0 from closed loop
        poly_proj = poly_t.projection(projs(p));
        
        poly_proj_vert = poly_proj.V';
        if length(poly_proj_vert) == 1
            % projection collapsed into a single point, duplicate for plot reasons
            disp("dim:" + length(poly_proj_vert))
            vert_projs{y,p} = repmat(poly_proj_vert,1,2);
        elseif isempty(poly_proj_vert)
            % MPT3 unable to find an interior due to numerical error. Warns
            % the user and use a backup method to find a feasible point.
            warning("dim:" + length(poly_proj_vert));
            poly_point = (poly_t.A(1:nx,:) \  poly_t.b(1:nx,:));
            vert_projs{y,p} = repmat(poly_point(projs(p)),1,2);
        else
            vert_projs{y,p} = poly_proj_vert;
        end
    end
end

end

function longDyn = defLongDyn(Ts,vx_interval,N_mpc,r_change)
% we define an acceleration phase (starting from t = r_change), going from
% vx_min to vx_max.
longDyn = zeros(2,N_mpc);
longDyn(2,1) = min(vx_interval);
for t = 1:N_mpc
    if t >= r_change
        vx_ref = max(vx_interval);
    else
        vx_ref = min(vx_interval);
    end
    if t < N_mpc % just avoid last time step
        vx_err = longDyn(2,t) - vx_ref;
        u_sym = -(0.4)*vx_err - sign(vx_err); % 6500 N / 2000 kg
        longDyn(:,t+1) = (eye(2)+Ts*[0,1;0,0])*longDyn(:,t) + Ts*[0;1]*u_sym;
        if abs(vx_err) < 0.01
            longDyn(2,t+1) = vx_ref;
        end
    end
end
longDyn = struct("pos",longDyn(1,:),"vx",longDyn(2,:));
end

function [A_convh_curr, B_convh_curr] = updateDynModel(A_is, B_is, P_hat, vx_curr)
% we are monotonically reducing the vx range. For this reason, P_hat is
% changing as if we add a new tighter constraint, replacing the previous one

% NOTE: here we don't collapse the set to a single point, we leave 5km/h range
P_hat_curr = Polyhedron([P_hat.A;-1,0],[P_hat.b;-min(vx_curr,max(P_hat.V(:,1))-5/3.6 )]);
P_hat_V = num2cell(P_hat_curr.V',1); np = length(P_hat_V);

A_convh_curr = cell(np,1); B_convh_curr = cell(np,1);
for v=1:np
    p = [P_hat_V{v};1];
    Atmp = p(1)*A_is{1}+p(2)*A_is{2}+p(3)*A_is{3};
    Btmp = p(1)*B_is{1}+p(2)*B_is{2}+p(3)*B_is{3};
    A_convh_curr{v} = Atmp;
    B_convh_curr{v} = Btmp;
end
end