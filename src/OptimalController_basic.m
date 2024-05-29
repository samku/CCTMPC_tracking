classdef OptimalController_basic < handle
    % OCP as defined in Equation (19) from the paper. The main difference
    % from the OptimalController class is the use of "optimize" function 
    % instead of "optimizer" from YALMIP. This leads to a simpler but less 
    % efficient code, as DynSystem is an handle class whose system matrices
    % are updated across all the instances.
    
    properties (SetAccess = private)
        sys % dynamical system
        N % prediction horizon length
        ccPoly % Configuration-Constrained Polytope
        sdp_opts % YALMIP optimizer options (solver, verbosity...)
        costFunMan % support class for cost function definition
        gamma % user specified discount factor
        
        feasRegion % Polytope, (approximated) feasible region for the OCP
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
        function obj = OptimalController_basic(sys, ccPoly, costFunMan, N, gamma, sdp_opts)
            obj.sys = sys;
            obj.ccPoly = ccPoly;
            obj.costFunMan = costFunMan;
            obj.N = N;
            obj.gamma = gamma;
            obj.sdp_opts = sdp_opts;
            
            obj.feasRegion = nan; % computation (slow) delayed, as only used in plots
        end
        
        function feas_x0 = findFeasibleX0(obj, x0_des, costFun)
            % find the closest feasible initial state to the desired one
            
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            
            % OCP_optimizer YALMIP definitions
            x0=sdpvar(obj.sys.nx, 1,'full');
            y=sdpvar(obj.ccPoly.m, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');
            ys=sdpvar(obj.ccPoly.m, 1,'full');
            us=sdpvar(obj.sys.nu, obj.ccPoly.v, 1,'full');
            
            switch class(costFun)
                case 'double'
                    cost = weightSq2Norm(x0 - x0_des, costFun);
                case 'function_handle'
                    cost = costFun(x0, x0_des);
                otherwise
                    error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(Stage_cost))
            end
            
            % % Constraints (look Equation (19) from paper)
            constr = [];
            
            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];
            
            % Stage + terminal constraints
            for k=1:obj.N-1
                constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1) )];
            end
            constr = [constr; obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.gamma*y(:,obj.N)+(1-obj.gamma)*ys)];
            
            % RCI constraint
            constr = [constr; obj.constrSetS(ys, us, ys)];
            
            % solve the OCP
            optimize(constr,cost,sdpsettings(obj.sdp_opts{:}));
            feas_x0 = value(x0);
        end
        
        
        function ocpSol = OCP_optimizer(obj, x0, r)
            % compute the main QP solution
            y=sdpvar(obj.ccPoly.m, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');
            ys=sdpvar(obj.ccPoly.m, 1,'full');
            us=sdpvar(obj.sys.nu, obj.ccPoly.v, 1,'full');
            
            % % Cost Function
            theta_s = sdpvar(obj.sys.nu,1);
            cost = obj.costFunMan.cost_RCI(ys,us,r,theta_s);
            for t=1:obj.N-1
                cost = cost + obj.costFunMan.cost_L_k(y(:,t), u(:,:,t), ys, us);
            end
            cost = cost + obj.costFunMan.cost_L_N(y(:,obj.N), u(:,:,obj.N), ys, us);
            
            % % Constraints (look Equation (19) from paper)
            constr = [];
            
            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];
            
            % Stage + terminal constraints
            for k=1:obj.N-1
                constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1) )];
            end
            constr = [constr; obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.gamma*y(:,obj.N)+(1-obj.gamma)*ys)];
            
            % RCI constraint
            constr = [constr; obj.constrSetS(ys, us, ys)];
            
            % solve the OCP
            optimize(constr,cost,sdpsettings(obj.sdp_opts{:}));
            
            % return optimization problem variables and cost
            ocpSol = {value(y),value(u),value(ys),value(us),value(cost)};
        end
        
        
        function rciSol = RCI_optimizer(obj, r)
            % compute the RCI set solution
            ys = sdpvar(obj.ccPoly.m, 1,'full');
            us = sdpvar(obj.sys.nu, obj.ccPoly.v,'full');
            
            % % Cost function
            theta_s = sdpvar(obj.sys.nu, 1,'full');
            cost = obj.costFunMan.cost_RCI(ys, us, r, theta_s);
            
            % RCI constraint
            constr = obj.constrSetS(ys,us,ys);
            
            % solve the OCP
            optimize(constr,cost,sdpsettings(obj.sdp_opts{:}));
            
            % return optimization problem variables and cost
            rciSol = {value(ys),value(us),value(cost)};
        end
        
        
        function feasRegion = computeFeasRegion(obj, n_facets)
            % compute a new feasible region given a template of #facets
            % complexity.
            H_feas = getTemplateMatrix(obj.sys.nx, n_facets);
            h_feas = zeros(size(H_feas,1),1);
            
            x0 = sdpvar(obj.sys.nx,1,'full');
            y=sdpvar(obj.ccPoly.m, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');
            ys=sdpvar(obj.ccPoly.m, 1,'full');
            us=sdpvar(obj.sys.nu, obj.ccPoly.v, 1,'full');
            
            for i=1:size(H_feas,1)                
                % cost function
                cost = -H_feas(i,:)*x0; % maximize the direction
                
                % % Constraints (look Equation (19) from paper)
                constr = [];
                
                % Initial condition
                constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];
                % Stage + terminal constraints
                for k=1:obj.N-1
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1) )];
                end
                constr = [constr; obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.gamma*y(:,obj.N)+(1-obj.gamma)*ys)];
                % RCI constraint
                constr = [constr; obj.constrSetS(ys, us, ys)];
                
                % solve the OCP
                optimize(constr,cost,sdpsettings(obj.sdp_opts{:}));
                
                % return the height on the specific direction
                h_feas(i) = -value(cost);
            end
                        
            feasRegion = Polyhedron(H_feas, h_feas).minHRep().minVRep();
        end
        
    end
    
    methods (Access = private)
        
        function constr = constrSetS(obj,y,u,yp)
            % define the Configuration Constrained RFIT set.
            % NOTE: if the system convex hulls are updated, the constraints
            % do it accordingly
            constr = [];
            for j=1:obj.ccPoly.v
                for i=1:obj.sys.np % #models
                    constr = [constr;
                        obj.ccPoly.F*(obj.sys.A_convh{i}*obj.ccPoly.Vi_s{j}*y + obj.sys.B_convh{i}*u(:,j) )+ obj.ccPoly.d <= yp];
                end
                constr = [constr;   obj.sys.X.A*obj.ccPoly.Vi_s{j}*y <= obj.sys.X.b;    obj.sys.U.A*u(:,j) <= obj.sys.U.b];
                
            end
            constr = [constr; obj.ccPoly.E*y <= 0];
        end
        
    end
    
end