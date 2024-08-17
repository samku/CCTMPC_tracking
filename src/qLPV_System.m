classdef qLPV_System < handle
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

    properties (Access = private)
        A_curr_han, B_curr_han  % function handle to get current system realization
    end

    properties(Dependent)
        nx, nu, ny  % state input and output dimensions
        np          % #vertices of convex hull
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
        function obj = qLPV_System(A_convh, B_convh, C, D, X, U, W_dist, A_curr_han, B_curr_han)
            obj.A_convh = A_convh;
            obj.B_convh = B_convh;
            obj.C = C;
            obj.D = D;
            obj.X = X;
            obj.U = U;
            obj.W_dist = W_dist;
            obj.A_curr_han = A_curr_han;
            obj.B_curr_han = B_curr_han;
        end

        function updateSysMatrices(obj, A_convh_new, B_convh_new)
            % used for online model refinement
            obj.A_convh = A_convh_new;
            obj.B_convh = B_convh_new;
        end

        function x_next = step_nominal(obj, x, u, varargin)
            % varargin handles external parameter vector for an LPV system
            x_next = obj.A_curr(x,varargin{:})*x ...
                + obj.B_curr(x,varargin{:})*u;
        end

        function A_curr = A_curr(obj,varargin)
            % varargin: {state, ext_param_vector}
            A_curr = obj.A_curr_han(varargin{:});
        end

        function B_curr = B_curr(obj,varargin)
            % varargin: {state, ext_param_vector}
            B_curr = obj.B_curr_han(varargin{:});
        end

    end

end