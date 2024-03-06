addpath(genpath('../src'))

% ---------------
% Note: we require that those softwares are already installed on the system
addpath(genpath('~/workSW/gurobi1100/linux64/matlab/'))
run('~/Documenti/CCTMPC/tbxmanager/startup.m')
addpath('~/workSW/casadi-3.6.4-linux64-matlab2018b/')
import casadi.*
% ---------------


% System definition
A=(1/5)*[5 1; -1 4];    nx = size(A,2);
B=(1/5)*[0;1];          nu = size(B,2);
Bw=(1/5)*[1 0; 0 1];    nw = size(Bw,2);

HX = [eye(nx); -eye(nx)]; hX = [6.8;10; 10;5];  mX = size(HX,1);
HU = [eye(nu); -eye(nu)]; hU = [10; 10];        mU = size(HU,1);
HW = [eye(nx); -eye(nx)]; hW = [0.5;2; 0.5;2];  mW = size(HW,1);

% output function
C = [1, 0]; ny = size(C,1);

% LTI system
A_convh = {}; B_convh = {};
A_convh{end+1} = A; B_convh{end+1} = B;

X = Polyhedron(HX, hX);
U = Polyhedron(HU, hU);
W_dist = Bw*Polyhedron(HW, hW);

dynSystem = DynSystem(A_convh, B_convh, X, U, W_dist);

% get Configuration-Constrained polytope
C_tilde = getConfConstrTemplate(nx,4);
ccpoly = CCPolytope(dynSystem, C_tilde, C_tilde);


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