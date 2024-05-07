classdef CCTMPC < handle
    %CCTMPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        var_convh % Variable system matrices flag
        
        r_prev % reference at previous timestep. Needed to trigger RCI computation
        ocpSol % cell array, solution from the main QP problem.
        rciSol % cell array, solution from RCI
        lambdaSol % solution of convex combination of the vertex control law
    end
    
    properties (Access = private)
        optController % YALMIP optimizer, Optimal Controller
        mpcLaw % YALMIP optimizer, output the input in state (and not parameter) space
    end
    
    properties(Dependent)
        Lyapunov_cost
        feasRegion % Polytope, define the (approximated) feasible region for the MPC scheme
    end
    
    methods % GETTER methods
        function Lyap_cost = get.Lyapunov_cost(obj)
            Lyap_cost = obj.ocpSol{end} - obj.rciSol{end};
        end
        
        function feasReg = get.feasRegion(obj)
            feasReg = obj.optController.feasRegion;
        end
    end
    
    
    methods (Access = public)
        function obj = CCTMPC(sys, ccPoly, CFM, N, gamma, sdp_opts, var_convh, useBasicOCP)
            %CCTMPC Construct an instance of this class
            %   Detailed explanation goes here
            if useBasicOCP
                obj.optController = OptimalController_basic(sys, ccPoly, CFM, N, gamma, sdp_opts);
            else
                obj.optController = OptimalController(sys, ccPoly, CFM, N, gamma, sdp_opts, var_convh);
                
            end
            obj.mpcLaw = obj.initMPCLaw();
            
            obj.var_convh = var_convh; % Variable system matrices flag
            
            obj.r_prev = nan;
            obj.ocpSol = cell(5,1); % {y,u,ys,us,cost} fields
            obj.rciSol = cell(3,1); % {ys,us,cost} fields
        end
        
        
        function uMPC = solve(obj, x_curr, r_curr, varargin)
            % solve OCP (in parameter space)
            switch class(obj.optController)
                case 'OptimalController'
                    obj.solve_fast(x_curr, r_curr, varargin{:});
                case 'OptimalController_basic'
                    obj.solve_basic(x_curr, r_curr, varargin{:});
                otherwise
                    error("Unsupported OCP class: %s", class(obj.optController))
            end
            % compute the optimal input (in state space) and lambda multipliers
            mpcSol = obj.mpcLaw(x_curr, obj.ocpSol{1}(:,1), obj.ocpSol{2}(:,:,1));
            [uMPC, obj.lambdaSol] = mpcSol{:};
            if isnan(uMPC)
                if sum(abs(obj.findFeasibleX0(x_curr,eye(obj.optController.sys.nx)) - x_curr)) > 1e-5
                    error("The current state does not lie inside the feasible"+...
                        " region of the control scheme. Consider changing initial condition.")
                end
                error("uMPC returned NaN. Consider enabling debug mode.")
            end
        end
        

        function feasRegion = computeFeasRegion(obj, n_facets, varargin)
            feasRegion = obj.optController.computeFeasRegion(n_facets, varargin{:});
        end
        
        
        function feasX0 = findFeasibleX0(obj,x_curr,varargin)
            % interface method to the underlying OCP
            switch length(varargin)
                case 0
                    feasX0 = obj.optController.findFeasibleX0(x_curr,eye(obj.optController.sys.nx));
                case 1
                    feasX0 = obj.optController.findFeasibleX0(x_curr,varargin{:});
                otherwise
                    error("%d arguments provided (max 1 expected).",length(varargin));
            end
        end
    end
    
    methods (Access = private)
        
        function solve_fast(obj, x_curr, r_curr, varargin)
            % two possible varargin cases:
            % length(varargin)==0: for fixed systems, RCI computed only if
            %       r_curr is different from previous one
            % length(varargin)==2: system dynamics (A_convh,B_convh) are
            %       updated; RCI always computed.
            switch length(varargin)
                case 0
                    if ~isequal(obj.r_prev, r_curr) % different reference
                        obj.r_prev = r_curr;
                        % update RCI set
                        if obj.var_convh
                            obj.rciSol = obj.optController.RCI_optimizer(r_curr,obj.optController.sys.A_convh,obj.optController.sys.B_convh);
                        else
                            obj.rciSol = obj.optController.RCI_optimizer(r_curr);
                        end
                    end
                    % solve OCP (solution in parameter space)
                    if obj.var_convh
                        obj.ocpSol = obj.optController.OCP_optimizer(x_curr, r_curr,obj.optController.sys.A_convh,obj.optController.sys.B_convh);          
                    else
                        obj.ocpSol = obj.optController.OCP_optimizer(x_curr, r_curr);
                    end
                case 2 % (A_convh,B_convh)
                    assert(obj.var_convh == true,...
                        "Dynamic system matrices have been initialized as fixed, cannot be modified online!");
                    obj.optController.sys.updateSysMatrices(varargin{:});
                    obj.r_prev = r_curr;
                    % update RCI set
                    obj.rciSol = obj.optController.RCI_optimizer(r_curr, varargin{:});
                    % solve OCP (solution in parameter space)
                    obj.ocpSol = obj.optController.OCP_optimizer(x_curr, r_curr,varargin{:});
                    
                otherwise % nargin-1 as obj is counted as argument
                    error("Wrong number of arguments: %d. Only 2 and 4 supported.", nargin-1);
            end
        end
        
        
        function solve_basic(obj, x_curr, r_curr, varargin)
            % two possible varargin cases:
            % length(varargin)==0: for fixed systems, RCI computed only if
            %       r_curr is different from previous one
            % length(varargin)==2: system dynamics (A_convh,B_convh) are
            %       updated; RCI always computed.
            
            % nargin-1 as obj is counted as argument
            assert(any(length(varargin)==[0,2]),...
                "Wrong number of arguments: %d. Only 2 and 4 supported.", nargin-1)
            
            % update system matrices
            if length(varargin) == 2
                obj.optController.sys.updateSysMatrices(varargin{:});
            end
            
            % if we have different reference or updating matrices
            if ~isequal(obj.r_prev, r_curr) || length(varargin) == 2
                obj.r_prev = r_curr;
                % update RCI set
                obj.rciSol = obj.optController.RCI_optimizer(r_curr);
            end
            
            % solve OCP (solution in parameter space)
            obj.ocpSol = obj.optController.OCP_optimizer(x_curr, r_curr);
        end
        
        
        function mpcLaw = initMPCLaw(obj)
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            
            m = obj.optController.ccPoly.m;
            m_bar = obj.optController.ccPoly.m_bar;
            nu = obj.optController.sys.nu;
            nx = obj.optController.sys.nx;
            Vi_s = obj.optController.ccPoly.Vi_s;
            
            y0Opt = sdpvar(m, 1,'full');
            u0Opt = sdpvar(nu, m_bar,'full');
            x0 = sdpvar(nx, 1,'full');
            
            lambda = sdpvar(m_bar, 1,'full');
            
            constr = [];
            constr = [constr; sum(lambda)==1; lambda >= 0];
            x_eq_sum = 0;
            for j=1:m_bar
                x_eq_sum = x_eq_sum + lambda(j)*Vi_s{j}*y0Opt;
            end
            constr = [constr; x0 == x_eq_sum];
            
            % cost = weightSq2Norm(lambda, eye(size(lambda,1))); % squared 2-norm
            cost = nnz(lambda);
            
            mpcLaw = optimizer(constr, cost, ...
                sdpsettings(obj.optController.sdp_opts{:}), ...
                {x0, y0Opt, u0Opt}, ... % input parameters
                {u0Opt*lambda, lambda} ... % output parameters
                ); % (nu x m_bar) * (m_bar x 1) -> (nu x 1)
        end
        
    end
    
end