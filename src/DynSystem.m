classdef DynSystem < handle
    %DYNSYS Summary of this class goes here
    %   Detailed explanation goes here
    
    % The system is expressed as x+ = Ax+Bu+w, where (A,B) \in \Delta,
    % convex hull of different realizations.
    
    properties (SetAccess = private)
        A_convh % Convex hull of state matrices
        B_convh % Convex hull of input matrices
        C, D    % state and input output matrices
        X       % State constraints
        U       % Input constraints 
        W_dist  % Disturbance constraints
    end
    
    properties(Dependent)
        nx, nu, ny  % state input and output dimensions
        np          % #vertices of convex hull
        
        % TODO: get actual realization!
    end
    
    methods % GETTER methods
        function nx = get.nx(obj)
            nx = size(obj.A_convh{1},2);
        end
        function nu = get.nu(obj)
            nu = size(obj.B_convh{1},2);
        end
        function ny = get.ny(obj)
            ny = size(obj.C,1);
        end
        function np = get.np(obj)
            np = length(obj.A_convh);
        end
    end
    
    methods (Access = public)
        function obj = DynSystem(A_convh, B_convh, C, D, X, U, W_dist)
            %DYNSYS Construct an instance of this class
            %   Detailed explanation goes here
            obj.A_convh = A_convh;
            obj.B_convh = B_convh;
            obj.C = C; 
            obj.D = D;
            obj.X = X;
            obj.U = U;
            obj.W_dist = W_dist;
        end
        
        function updateSysMatrices(obj, A_convh_new, B_convh_new)
        % used for online model refinement
            obj.A_convh = A_convh_new;
            obj.B_convh = B_convh_new;
        end
        
        function x_next = step(obj, x0, u0)
            % TODO: define A_curr and B_curr! maybe passing a function
            % handle to define how to compute them?
            x_next = A_curr*x0+ B_curr*u0 + obj.W_dist.randomPoint();
        end
        
    end
end

