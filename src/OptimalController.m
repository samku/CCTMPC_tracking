classdef OptimalController < handle
    % The optimal controller defind with YALMIP. preferred solver is gurobi
    % we also define a flag variable to disable matrix update if not
    % required and speeding up computation (2-3x improvement)
    
    properties (SetAccess = private)
        sys % dynamical system
        N % prediction horizon length
        ccPoly % Configuration constrained poly
        sdp_opts % options related to SDPVar (solver, verbosity...)
        costFunMan % support class for cost function definition
        OCP_optimizer  % main QP problem.
        RCI_optimizer % QP for solving RCI set only. Used mainly for plotting reasons
        gamma % user specified discount factor
        
        feasRegion % Polytope, define the (approximated) feasible region for the OCP scheme
    end
    
    properties (Access = private)
        directionalFeasRegion % YALMIP optimizer, compute the farthest feasible
        % initial condition with respect to a specific facet.
    end
    
    methods % GETTER methods
        function feasReg = get.feasRegion(obj)
            if class(obj.feasRegion) == "Polyhedron"
                feasReg = obj.feasRegion;
            else
                % Hard-coded #facets, often enough to get a good approximation
                obj.feasRegion = obj.computeFeasRegion(100);
                feasReg = obj.feasRegion;
            end
        end
    end
    
    methods (Access = public)
        function obj = OptimalController(sys, ccPoly, costFunMan, N, gamma, sdp_opts, var_convh)
            %OPTIMALCONTROLLER Construct an instance of this class
            %   Detailed explanation goes here
            
            % TODO: the var_convh is not sufficiently robust yet. Edge
            % cases must still addressed.
            
            obj.sys = sys;
            obj.N = N;
            obj.ccPoly = ccPoly;
            obj.sdp_opts = sdp_opts;
            obj.costFunMan = costFunMan;
            obj.gamma = gamma;
            obj.OCP_optimizer = obj.initOCPOptimizer(var_convh);
            obj.RCI_optimizer = obj.initRCIOptimizer(var_convh);
            
            obj.directionalFeasRegion = obj.initDirFeasRegion(var_convh);
            obj.feasRegion = nan; % computation delayed (slow), used only for plots.
        end
        
        
        function feas_x0 = findFeasibleX0(obj, x0_des, costFun)
            
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            
            % OCP_optimizer yalmip definitions
            x0=sdpvar(obj.sys.nx, 1,'full');
            y=sdpvar(obj.ccPoly.m, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.m_bar, obj.N,'full');
            ys=sdpvar(obj.ccPoly.m, 1,'full');
            us=sdpvar(obj.sys.nu, obj.ccPoly.m_bar, 1,'full');
            
            switch class(costFun)
                case 'double'
                    cost = weightSq2Norm(x0 - x0_des, costFun);
                case 'function_handle'
                    cost = costFun(x0, x0_des);
                otherwise
                    error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(Stage_cost))
            end
            
            % % Constraints (look Equation 17 from paper)
            constr = [];
            
            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];
            
            % Stage + terminal constraints
            for k=1:obj.N-1
                constr = [constr; 
                    obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1),obj.sys.A_convh,obj.sys.B_convh )];
            end
            constr = [constr; 
                obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.gamma*y(:,obj.N)+(1-obj.gamma)*ys,...
                obj.sys.A_convh,obj.sys.B_convh)];
            
            % RCI constraint
            constr = [constr; obj.constrSetS(ys, us, ys,obj.sys.A_convh,obj.sys.B_convh)];
            
            % solve the OCP
            optimize(constr,cost,sdpsettings(obj.sdp_opts{:}));
            feas_x0 = value(x0);
        end
        
        
        function feasRegion = computeFeasRegion(obj, n_facets, varargin)
            % compute a new feasible region for a template of #facets
            % complexity. varargin used to pass new (A_convh, B_convh)
            H_feas = getTemplateMatrix(obj.sys.nx, n_facets);
            h_feas = zeros(size(H_feas,1),1);
            
            for i=1:size(H_feas,1)
                switch length(varargin)
                    case 0
                        x0_max = obj.directionalFeasRegion(H_feas(i,:)');
                    case 2
                        x0_max = obj.directionalFeasRegion(H_feas(i,:)',varargin{:});
                    otherwise % nargin-1 as obj is counted as argument
                        error("Wrong number of arguments: %d. Only 2 and 4 supported.", nargin-1);
                end
                h_feas(i) = H_feas(i,:)*x0_max;
            end
            
            feasRegion = Polyhedron(H_feas, h_feas).minHRep().minVRep();
        end
    end
    
    
    methods (Access = private)        
        function ocpOptim = initOCPOptimizer(obj, var_convh)
            y=sdpvar(obj.ccPoly.m, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.m_bar, obj.N,'full');
            x0=sdpvar(obj.sys.nx, 1,'full');
            ys=sdpvar(obj.ccPoly.m, 1,'full');
            us=sdpvar(obj.sys.nu, obj.ccPoly.m_bar, 1,'full');
            r=sdpvar(obj.sys.ny ,1,'full');
            
            if var_convh
                A_sdp = {};
                for i = 1:obj.sys.np
                    A_sdp{i} = sdpvar(obj.sys.nx,obj.sys.nx,'full');
                end
                B_sdp = {};
                for i = 1:obj.sys.np
                    B_sdp{i} = sdpvar(obj.sys.nx,obj.sys.nu,'full');
                end
            end
            
            % % Cost Function
            theta_s = sdpvar(obj.sys.nu,1);
            cost = obj.costFunMan.cost_RCI(ys,us,r,theta_s);
            for t=1:obj.N-1
                cost = cost + obj.costFunMan.cost_L_k(y(:,t), u(:,:,t), ys, us);
            end
            cost = cost + obj.costFunMan.cost_L_N(y(:,obj.N), u(:,:,obj.N), ys, us);
            
            % % Constraints (look Equation 17 from paper)
            constr = [];
            
            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];
            
            % Stage + terminal constraints
            for k=1:obj.N-1
                if var_convh
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1), A_sdp,B_sdp )];
                else
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1), obj.sys.A_convh,obj.sys.B_convh )];
                end
            end
            if var_convh
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.gamma*y(:,obj.N)+(1-obj.gamma)*ys, A_sdp,B_sdp)];
            else
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.gamma*y(:,obj.N)+(1-obj.gamma)*ys,  obj.sys.A_convh,obj.sys.B_convh)];
            end
            
            % RCI constraint
            if var_convh
                constr = [constr; obj.constrSetS(ys, us, ys,  A_sdp,B_sdp)];
            else
                constr = [constr; obj.constrSetS(ys, us, ys,  obj.sys.A_convh,obj.sys.B_convh)];
            end
            
            if var_convh
                ocpOptim = optimizer(constr,cost,...
                    sdpsettings(obj.sdp_opts{:}),...
                    {x0,r,A_sdp{:},B_sdp{:}},...% input parameters
                    {y,u,ys,us,cost});          % output parameters
            else
                ocpOptim = optimizer(constr,cost,...
                    sdpsettings(obj.sdp_opts{:}),...
                    {x0,r},...          % input parameters
                    {y,u,ys,us,cost});  % output parameters
            end
        end
        
        
        function rciOptim = initRCIOptimizer(obj, var_convh)
            ys = sdpvar(obj.ccPoly.m, 1,'full');
            us = sdpvar(obj.sys.nu, obj.ccPoly.m_bar,'full');
            r = sdpvar(obj.sys.ny, 1,'full');
            
            if var_convh
                A_sdp = {};
                for i = 1:obj.sys.np
                    A_sdp{i} = sdpvar(obj.sys.nx,obj.sys.nx,'full');
                end
                B_sdp = {};
                for i = 1:obj.sys.np
                    B_sdp{i} = sdpvar(obj.sys.nx,obj.sys.nu,'full');
                end
            end
            
            % % Cost function
            theta_s = sdpvar(obj.sys.nu, 1,'full');
            cost = obj.costFunMan.cost_RCI(ys, us, r, theta_s);
            
            % RCI constraint
            if var_convh
                constr = obj.constrSetS(ys,us,ys, A_sdp,B_sdp);
            else
                constr = obj.constrSetS(ys,us,ys, obj.sys.A_convh, obj.sys.B_convh);
            end
            
            if var_convh
                rciOptim = optimizer(constr, cost,...
                    sdpsettings(obj.sdp_opts{:}),...
                    {r,A_sdp{:},B_sdp{:}},... % input parameters
                    {ys,us,cost});          % output parameters
            else
                rciOptim = optimizer(constr, cost,...
                    sdpsettings(obj.sdp_opts{:}),...
                    {r},...         % input parameters
                    {ys,us,cost});  % output parameters
            end
        end
        
        
        function dirFeasRegion = initDirFeasRegion(obj, var_convh)
            c_vec = sdpvar(obj.sys.nx,1,'full');
            x0 = sdpvar(obj.sys.nx,1,'full');
            
            % cost function
            cost = -c_vec'*x0; % maximize the direction
            
            y=sdpvar(obj.ccPoly.m, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.m_bar, obj.N,'full');
            ys=sdpvar(obj.ccPoly.m, 1,'full');
            us=sdpvar(obj.sys.nu, obj.ccPoly.m_bar, 1,'full');
            
            if var_convh
                A_sdp = {};
                for i = 1:obj.sys.np
                    A_sdp{i} = sdpvar(obj.sys.nx, obj.sys.nx,'full');
                end
                B_sdp = {};
                for i = 1:obj.sys.np
                    B_sdp{i} = sdpvar(obj.sys.nx, obj.sys.nu,'full');
                end
            end
            
            % % Constraints (look Equation 17 from paper)
            constr = [];
            
            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];
            
            % Stage + terminal constraints
            for k=1:obj.N-1
                if var_convh
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1), A_sdp,B_sdp )];
                else
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1), obj.sys.A_convh,obj.sys.B_convh )];
                end
            end
            if var_convh
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.gamma*y(:,obj.N)+(1-obj.gamma)*ys, A_sdp,B_sdp)];
            else
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.gamma*y(:,obj.N)+(1-obj.gamma)*ys,  obj.sys.A_convh,obj.sys.B_convh)];
            end
            
            % RCI constraint
            if var_convh
                constr = [constr; obj.constrSetS(ys, us, ys,  A_sdp,B_sdp)];
            else
                constr = [constr; obj.constrSetS(ys, us, ys,  obj.sys.A_convh,obj.sys.B_convh)];
            end
            
            if var_convh
                dirFeasRegion = optimizer(constr, cost,...
                    sdpsettings(obj.sdp_opts{:}),...
                    {c_vec,A_sdp{:},B_sdp{:}},... % input parameters
                    {x0});      % output parameters
            else
                dirFeasRegion = optimizer(constr, cost,...
                    sdpsettings(obj.sdp_opts{:}),...
                    {c_vec},... % input parameters
                    {x0});      % output parameters
            end
        end
        
        
        function constr = constrSetS(obj,y,u,yp,A_convh,B_convh)
            % define the Configuration Constrained RFIT set.
            constr = [];
            for j=1:obj.ccPoly.m_bar
                for i=1:obj.sys.np % #models
                    constr = [constr;
                        obj.ccPoly.F*(A_convh{i}*obj.ccPoly.Vi_s{j}*y + B_convh{i}*u(:,j) )+ obj.ccPoly.d <= yp];
                end
                constr = [constr;   obj.sys.X.A*obj.ccPoly.Vi_s{j}*y <= obj.sys.X.b;    obj.sys.U.A*u(:,j) <= obj.sys.U.b];
                
            end
            constr = [constr; obj.ccPoly.E*y <= 0];
        end
        
    end
    
end