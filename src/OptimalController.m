classdef OptimalController < handle
    %OPTIMALCONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        sys % dynamical system
        N
        ccPoly
        solver
        costFunMan
        optContr
        optRCI
    end
    
    methods (Access = public)
        function obj = OptimalController(dynSys, ccPoly, costFunMan, N, gamma, solver_opts)
            %OPTIMALCONTROLLER Construct an instance of this class
            %   Detailed explanation goes here
            
            % again, varargin to not specify solver eventually.
            % TODO: the cost must be defined inside this class? as it
            % depends on the number of steps to compute the more complex
            % one. Or instead building another class, and share N for both
            % OptContr and Cost classes on the script example?
            
            obj.sys = dynSys;
            obj.N = N;
            obj.ccPoly = ccPoly;
            obj.solver = solver_opts;
            
            obj.optContr = obj.buildOptController(costFunMan, N, gamma, solver_opts);
            obj.optRCI = obj.buildOptRCI(costFunMan, solver_opts);
            
        end
        
        
        
        
    end
    
    methods (Access = private)
        
        function optController = buildOptController(obj, CFM, N, gamma, solver_opts)
            y=sdpvar(obj.ccPoly.m, N,'full');
            u=sdpvar(nu, obj.ccPoly.m_bar, N,'full');
            x0=sdpvar(nx, 1,'full');
            ys=sdpvar(obj.ccPoly.m, 1,'full');
            us=sdpvar(nu, obj.ccPoly.m_bar, 1,'full');
            r=sdpvar(ny ,1,'full');
            
            % % Cost Function
                theta_s = sdpvar(nu,1);
                cost = CFM.cost_RCI(ys,us,r,theta_s);
                for t=1:N-1
                    cost = cost + CFM.cost_L_k(y(:,t), u(:,:,t), ys, us);
                end
                cost = cost + CFM.cost_L_N(y(:,N), u(:,:,N), ys, us);
            
            % % Constraints (look Problem 8 and Equation 17 from paper)
            constr = [];
            
            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];
            
            % Stage + terminal constraints
            for k = 1:N
                for j=1:obj.ccPoly.m_bar
                    for i=1:obj.sys.np % #models
                        if k < N
                            constr = [constr;
                                obj.ccPoly.F*(A_convh{i}*Vi_s{j}*y(:,k) + B_convh{i}*u(:,j,k)) + d <= y(:,k+1)];
                        else % last condition
                            constr = [constr;
                                obj.ccPoly.F*(A_convh{i}*Vi_s{j}*y(:,N) + B_convh{i}*u(:,j,N)) + d <= gamma*y(:,N) + (1-gamma)*ys];
                        end
                    end
                    constr = [constr;   U.A*u(:,j,k) <= U.b;    X.A*Vi_s{j}*y(:,k) <= X.b];
                end
                constr = [constr; obj.ccPoly.E*y(:,k) <= zeros(size(obj.ccPoly.E,1),1)];
            end
            
            % RCI constraint
            for j=1:obj.ccPoly.m_bar
                for i=1:obj.sys.np % #models
                    constr = [constr;
                        obj.ccPoly.F*(A_convh{i}*Vi_s{j}*ys + B_convh{i}*us(:,j)) + d <= ys];
                end
                constr = [constr;   U.A*us(:,j) <= U.b;    X.A*Vi_s{j}*ys <= X.b];
            end
            constr = [constr; obj.ccPoly.E*ys <= zeros(size(obj.ccPoly.E,1),1)];
            
            optController = optimizer(constr,cost,...
                sdpsettings(solver_opts{:}),...
                {x0,r},...                  % input parameters
                {y,u,cost,ys,us,theta_s});  % output parameters
            
        end
        
        function optRCI = buildOptRCI(obj, CFM, solver_opts)
            ys = sdpvar(obj.ccPoly.m, 1,'full');
            us = sdpvar(obj.sys.nu, obj.ccPoly.m_bar,'full');
            r = sdpvar(obj.sys.ny, 1,'full');
            
            % % Cost function
            theta_s = sdpvar(obj.sys.nu, 1,'full');
            cost = CFM.cost_RCI(ys, us, r, theta_s);
            
            % RCI constraint
            constr = [];
            for j=1:obj.ccPoly.m_bar
                for i=1:obj.sys.np % #models
                    constr = [constr;
                        obj.ccPoly.F*(obj.sys.A_convh{i}*obj.ccPoly.Vi_s{j}*ys + obj.sys.B_convh{i}*us(:,j)) + obj.ccPoly.d <= ys];
                end
                constr = [constr;   obj.sys.U.A*us(:,j) <= obj.sys.U.b;    obj.sys.X.A*obj.ccPoly.Vi_s{j}*ys <= obj.sys.X.b];
            end
            constr = [constr; obj.ccPoly.E*ys <= zeros(size(obj.ccPoly.E,1),1)];
            
            optRCI = optimizer(constr, cost,...
                sdpsettings(solver_opts{:}),...
                {z_r},...               % input parameters
                {ys,us,cost,theta_s});  % output parameters
        end
        
        
        function RCI_opt = computeRCI(obj,z_r)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            % TODO: add a flag such that only if specified we care about variable
            % state and input matrices!
            
            RCI_opt = 0;
        end
        
        
    end
end