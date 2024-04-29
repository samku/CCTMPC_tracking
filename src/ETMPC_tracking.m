classdef ETMPC_tracking < handle
    % ETMPC_TRACKING Summary of this class goes here
    %   Detailed explanation goes here
    % ETMPC implementation has been modified according to "RTMPC for
    % tracking" by Limon. In this way, we increase feasible region by
    % targeting an artificial output.
    
    properties (SetAccess = private)
        optContrSol % struct, solution from OCP
        
        feasRegion % Polytope, define the (approximated) feasible region for the MPC scheme
        
        sys % system dynamics and constraints
        ccPoly % vertex and parameters from the CCPolytope are used for comparison
        costMats % struct that aggregates all the relevant matrices for the system
        N % prediction horizon for the OCP
        lambda % user defined discount factor for the terminal set
        sdp_opts % options related to SDPVar (solver, verbosity...)
        
        y_mRPI
        L_RPI, L_HX, L_HU
        termSetAug
        K, M, P % affine feedback law, nullspace of steady states and terminal state cost matrix
        Q_RPI, P_RPI % cost matrices for the elastic term a_t
    end
    
    properties (Access = private)
        optController % YALMIP optimizer, Optimal Controller
        K_projector % YALMIP optimizer, for closed loop input constraint satisfaction
        directionalFeasRegion % YALMIP optimizer, compute the farthest feasible
        % initial condition with respect to a specific facet.
    end
    
    properties(Dependent)
        Lyapunov_cost
    end
    
    methods % GETTER methods
        function Lyap_cost = get.Lyapunov_cost(obj)
            Lyap_cost = obj.optContrSol.cost;
        end
        
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
        function obj = ETMPC_tracking(sys, ccPoly, costMatricesSS, N, lambda, sdp_opts)
            assert(sys.np == 1, "Cannot proceed. ETMPC from Limon is defined for only LTI systems.")
            
            obj.sys = sys;
            obj.ccPoly = ccPoly;
            obj.N = N;
            obj.lambda = lambda;
            obj.sdp_opts = sdp_opts;
            
            % struct{Qv_(y/u)p, Qv_(x/u), Qc_(y/u)p, Qc_(x/u), Qr}
            obj.costMats = obj.computeCostMats(costMatricesSS); %{Qv, Qc, Qr}
            
            % compute terminal set for [z;theta;a]. Note that some class
            % attributes are instantiated inside the method.
            obj.termSetAug = obj.computeTerminalSet();
            
            % instantiate YALMIP optimizers
            obj.optController = obj.initOptController();
            obj.K_projector = obj.initKProjector();
            
            obj.optContrSol = cell(5,1); % {x,u,theta,a,cost} fields
            
            obj.directionalFeasRegion = obj.initDirFeasRegion();
            obj.feasRegion = nan; % computation delayed (slow), used only for plots.
            
        end
        
        function uMPC = solve(obj, x_curr, r_curr)
            % solve OCP (solution in state space)
            obj.optContrSol = obj.optController(x_curr, r_curr);
            
            % compute closed loop solution: u_cl = v_0 + K_proj*(x_curr - z_0)
            uMPC = obj.optContrSol{2}(:,1) + obj.K_projector(obj.K)*(x_curr - obj.optContrSol{1}(:,1));
            if isnan(uMPC)
                if sum(abs(obj.findFeasibleX0(x_curr,eye(obj.sys.nx)) - x_curr)) > 1e-5
                    error("The current state does not lie inside the feasible"+...
                        " region of the control scheme. Consider changing initial condition.")
                end
                error("uMPC returned NaN. Consider enabling debug mode.")
            end
        end
        
        function feas_x0 = findFeasibleX0(obj, x0_des, costFun)
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            Mx = obj.M(1:obj.sys.nx,:); Mu = obj.M(obj.sys.nx+1:end,:);
            
            % optContr yalmip definitions
            theta = sdpvar(obj.sys.nu, 1,'full');
            x0 = sdpvar(obj.sys.nx, 1, 'full');
            x = sdpvar(obj.sys.nx, obj.N,'full');
            u = sdpvar(obj.sys.nu, obj.N-1,'full');
            a = sdpvar(obj.ccPoly.m, obj.N,'full');
            
            % define cost function
            switch class(costFun)
                case 'double'
                    cost = weightSq2Norm(x0 - x0_des, costFun);
                case 'function_handle'
                    cost = costFun(x0, x0_des);
                otherwise
                    error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(Stage_cost))
            end
            
            % initial condition & steady-state constraints
            constr = [obj.ccPoly.F*(x0-x(:,1)) <= a(:,1); % initial condition
                obj.sys.X.A*Mx*theta <= obj.sys.X.b;
                obj.sys.U.A*Mu*theta <= obj.sys.U.b];
            
            % stage cost and constraints
            for t = 1:obj.N-1
                constr = [constr;   % invariance constraint
                    obj.ccPoly.F*(A*x(:,t)+B*u(:,t))+obj.L_RPI*a(:,t)+obj.ccPoly.d <= obj.ccPoly.F*x(:,t+1)+a(:,t+1)];
                % state and input constraints
                constr=[constr; obj.sys.X.A*x(:,t)+obj.L_HX*a(:,t) <= obj.sys.X.b];
                constr=[constr; obj.sys.U.A*u(:,t)+obj.L_HU*a(:,t) <= obj.sys.U.b];
                % nonnegative elastic parameter
                constr=[constr; a(:,t) >= 0];
            end
            % terminal cost and constraints
            constr = [constr; obj.termSetAug.A*[x(:,obj.N); theta; a(:,obj.N)] <= obj.termSetAug.b];
            
            % solve the OCP
            optimize(constr,cost,sdpsettings(obj.sdp_opts{:}));
            feas_x0 = value(x0);
        end
        
        
        function feasRegion = computeFeasRegion(obj, n_facets)
            % compute a new feasible region for a template of #facets
            % complexity.
            H_stab = getTemplateMatrix(obj.sys.nx, n_facets);
            h_stab = zeros(size(H_stab,1),1);
            
            for i=1:size(H_stab,1)
                x0_max = obj.directionalFeasRegion(H_stab(i,:)');
                h_stab(i) = H_stab(i,:)*x0_max;
            end
            
            feasRegion = Polyhedron(H_stab, h_stab).minHRep().minVRep();
        end
    end
    
    
    methods (Access = private)
        function termSetAug = computeTerminalSet(obj)
            % frequently used constants throughout function
            A = obj.sys.A_curr(); B = obj.sys.B_curr();
            nx = obj.sys.nx; nu = obj.sys.nu;
            HX = obj.sys.X.A; hX = obj.sys.X.b;
            HU = obj.sys.U.A; hU = obj.sys.U.b;
            m = obj.ccPoly.m;
            
            % compute mRPI facet configuration (as Limon). LQR solutions
            % are assigned inside this method
            obj.y_mRPI = obj.get_mRPI_facets();
            
            % solve dual problem to get multipliers
            [obj.L_RPI, obj.L_HX, obj.L_HU] = obj.getMultipliers();
            
            % compute terminal cost for the elastic term [a]
            obj.Q_RPI = eye(obj.ccPoly.m);
            obj.P_RPI = obj.getRPITerminalCost();
            
            % compute the steady-state subspace
            obj.M = null([A-eye(nx) B]);
            Mx = obj.M(1:nx,:); Mu = obj.M(nx+1:end,:);
            
            
            % extend state vector and constraint matrices to [z;theta;a]
            A_aug = [A+B*obj.K, B*Mu-B*obj.K*Mx, zeros(nx, m);
                zeros(nu,nx), eye(nu),zeros(nu, m);
                zeros(m,nx)     zeros(m,nu)  obj.L_RPI];
            
            d_aug = [zeros(nx+nu,1); obj.ccPoly.d];
            
            H_aug = [HX, zeros(size(HX,1),nu), obj.L_HX;
                HU*obj.K,  HU*Mu-HU*obj.K*Mx, obj.L_HU;
                zeros(size(HX,1),nx), HX*Mx,  zeros(size(HX,1),m);
                zeros(size(HU,1),nx), HU*Mu,  zeros(size(HU,1),m);
                zeros(m,nx),  zeros(m,nu), -eye(m)];
            
            h_aug = [hX;
                hU;
                obj.lambda*hX;
                obj.lambda*hU;
                zeros(m,1)];
            
            % iterate constraints through system dynamics to overapproximate
            % the terminal set from above (volume of the set diminishes at
            % each iteration)
            H_RPI = H_aug; h_RPI = h_aug;
            t = 0;
            while true
                t = t+1; %disp(t)
                H_RPI = [H_RPI; H_aug*A_aug^t];
                % tighten constraints to get an invariant set
                rhs_reduce = zeros(size(h_aug));
                for jj=t-1:-1:0
                    rhs_reduce = rhs_reduce + H_aug*(A_aug^jj)*d_aug;
                end
                h_RPI = [h_RPI; h_aug - rhs_reduce];
                
                % Verify invariance H(Ax+d)<=h for all x:Hx<=h
                new_h_RPI = h_RPI;
                for i=1:size(h_RPI,1)
                    xx = cplexlp(-(H_RPI(i,:)*A_aug)',H_RPI,h_RPI);
                    new_h_RPI(i,1) = H_RPI(i,:)*(A_aug*xx+d_aug);
                end
                
                if max(new_h_RPI-h_RPI) <= 1e-6 % suff. small improvement
                    break
                end
            end
            % now (H_RPI,h_RPI) can be used for the terminal set.
            termSetAug = Polyhedron(H_RPI,h_RPI).minHRep();
        end
        
        
        function optContr = initOptController(obj)
            % frequently used constants throughout function
            A = obj.sys.A_curr(); B = obj.sys.B_curr(); C = obj.sys.C; D = obj.sys.D;
            Mx = obj.M(1:obj.sys.nx,:); Mu = obj.M(obj.sys.nx+1:end,:);
            
            % define OCP variables
            r = sdpvar(obj.sys.ny, 1,'full');
            theta = sdpvar(obj.sys.nu, 1,'full');
            x0 = sdpvar(obj.sys.nx, 1, 'full');
            x = sdpvar(obj.sys.nx, obj.N,'full');
            u = sdpvar(obj.sys.nu, obj.N-1,'full');
            a = sdpvar(obj.ccPoly.m, obj.N,'full');
            
            
            % initial condition & steady-state constraints
            constr = [obj.ccPoly.F*(x0-x(:,1)) <= a(:,1); % initial condition
                obj.sys.X.A*Mx*theta <= obj.sys.X.b;
                obj.sys.U.A*Mu*theta <= obj.sys.U.b];
            % tracking cost
            cost = (r - [C, D]*obj.M*theta)'*obj.costMats.Qr*(r - [C, D]*obj.M*theta);
            
            % stage cost and constraints
            for t = 1:obj.N-1
                constr = [constr;   % invariance constraint
                    obj.ccPoly.F*(A*x(:,t)+B*u(:,t))+obj.L_RPI*a(:,t)+obj.ccPoly.d <= obj.ccPoly.F*x(:,t+1)+a(:,t+1)];
                % state and input constraints
                constr=[constr; obj.sys.X.A*x(:,t)+obj.L_HX*a(:,t) <= obj.sys.X.b];
                constr=[constr; obj.sys.U.A*u(:,t)+obj.L_HU*a(:,t) <= obj.sys.U.b];
                % nonnegative elastic parameter
                constr=[constr; a(:,t) >= 0];
                
                cost = cost + (x(:,t)-Mx*theta)'*obj.costMats.Qc_x*(x(:,t)-Mx*theta)...
                    + (u(:,t)-Mu*theta)'*obj.costMats.Qc_u*(u(:,t)-Mu*theta) + ...
                    (a(:,t)-obj.y_mRPI)'*obj.Q_RPI*(a(:,t)-obj.y_mRPI);
            end
            % terminal cost and constraints
            constr = [constr; obj.termSetAug.A*[x(:,obj.N); theta; a(:,obj.N)] <= obj.termSetAug.b];
            cost = cost + (x(:,obj.N)-Mx*theta)'*obj.P*(x(:,obj.N)-Mx*theta) + ...
                (a(:,obj.N)-obj.y_mRPI)'*obj.P_RPI*(a(:,obj.N)-obj.y_mRPI);
            
            % instantiate Optimizer
            optContr = optimizer(constr,cost,sdpsettings(obj.sdp_opts{:}),...
                {x0,r},...      % input parameters
                {x,u,theta,a,cost});  % output parameters
        end
        
        
        function K_projector = initKProjector(obj)
            % ensure that the closed loop control law {u_cl = v0+K*(x-z0)}
            % fulfills input constraints.
            K_des = sdpvar(obj.sys.nu, obj.sys.nx,'full');
            K_proj= sdpvar(obj.sys.nu, obj.sys.nx,'full');
            constr = [obj.L_RPI*obj.ccPoly.F == obj.ccPoly.F*(obj.sys.A_curr()+obj.sys.B_curr()*K_proj);
                obj.L_HU*obj.ccPoly.F == obj.sys.U.A*K_proj];
            cost = norm(K_proj-K_des,1);
            K_projector = optimizer(constr,cost,sdpsettings(obj.sdp_opts{:}),...
                {K_des},... % input parameters
                {K_proj});  % output parameters
        end
        
        
        function dirFeasRegion = initDirFeasRegion(obj)
            % compute the directional feasibility region with respect to a
            % specific direction.
            A = obj.sys.A_curr(); B = obj.sys.B_curr();
            Mx = obj.M(1:obj.sys.nx,:); Mu = obj.M(obj.sys.nx+1:end,:);
            c_vec = sdpvar(obj.sys.nx,1,'full');
            
            theta = sdpvar(obj.sys.nu,1,'full');
            x0 = sdpvar(obj.sys.nx,1,'full');
            x = sdpvar(obj.sys.nx,obj.N,'full');
            u = sdpvar(obj.sys.nu,obj.N-1,'full');
            a = sdpvar(obj.ccPoly.m, obj.N,'full');
            
            % cost function
            cost = -c_vec'*x0; % maximize the direction
            
            % initial condition & steady-state constraints
            constr = [obj.ccPoly.F*(x0-x(:,1)) <= a(:,1); % initial condition
                obj.sys.X.A*Mx*theta <= obj.sys.X.b;
                obj.sys.U.A*Mu*theta <= obj.sys.U.b];
            
            % stage cost and constraints
            for t = 1:obj.N-1
                constr = [constr;   % invariance constraint
                    obj.ccPoly.F*(A*x(:,t)+B*u(:,t))+obj.L_RPI*a(:,t)+obj.ccPoly.d <= obj.ccPoly.F*x(:,t+1)+a(:,t+1)];
                % state and input constraints
                constr=[constr; obj.sys.X.A*x(:,t)+obj.L_HX*a(:,t) <= obj.sys.X.b];
                constr=[constr; obj.sys.U.A*u(:,t)+obj.L_HU*a(:,t) <= obj.sys.U.b];
                % nonnegative elastic parameter
                constr=[constr; a(:,t) >= 0];
            end
            % terminal constraints
            constr = [constr; obj.termSetAug.A*[x(:,obj.N); theta; a(:,obj.N)] <= obj.termSetAug.b];
            
            dirFeasRegion = optimizer(constr, cost, sdpsettings(obj.sdp_opts{:}),...
                {c_vec},... % input parameters
                {x0});      % output parameters
        end
        
        
        function y_mRPI = get_mRPI_facets(obj)
            % Limon mRPI set
            A = obj.sys.A_curr(); B = obj.sys.B_curr();
            
            [obj.K, obj.P] = dlqr(A, B, obj.costMats.Qc_x, obj.costMats.Qc_u);
            obj.K = -obj.K;
            
            xs = sdpvar(obj.sys.nx, obj.ccPoly.m,'full');
            c_vec = sdpvar(obj.ccPoly.m, 1, 'full');
            
            constr = [];
            for i=1:obj.ccPoly.m
                constr = [constr; c_vec(i) <= obj.ccPoly.F(i,:)*(A+B*obj.K)*xs(:,i)];
                constr = [constr; obj.ccPoly.F*xs(:,i) <= c_vec + obj.ccPoly.d];
            end
            cost = -sum(c_vec);
            optimize(constr, cost, sdpsettings(obj.sdp_opts{:}));
            
            y_mRPI = value(c_vec) + obj.ccPoly.d;
        end
        
        
        function [L_RPI, L_HX, L_HU] = getMultipliers(obj)
            % Solve the dual LP problem to get multipliers
            A = obj.sys.A_curr(); B = obj.sys.B_curr();
            
            L_RPI = zeros(obj.ccPoly.m, obj.ccPoly.m);
            L_HX = zeros(size(obj.sys.X.A,1), obj.ccPoly.m);
            L_HU = zeros(size(obj.sys.U.A,1), obj.ccPoly.m);
            
            for i=1:size(L_RPI,1)
                lpSol = cplexlp(obj.y_mRPI,...
                    -eye(obj.ccPoly.m), zeros(obj.ccPoly.m,1), ...
                    obj.ccPoly.F', (obj.ccPoly.F(i,:)*(A+B*obj.K))' );
                L_RPI(i,:) = lpSol';
            end
            
            for i=1:size(L_HX,1)
                lpSol = cplexlp(obj.y_mRPI,...
                    -eye(obj.ccPoly.m), zeros(obj.ccPoly.m,1), ...
                    obj.ccPoly.F', obj.sys.X.A(i,:)');
                L_HX(i,:) = lpSol';
            end
            for i=1:size(L_HU,1)
                lpSol = cplexlp(obj.y_mRPI,...
                    -eye(obj.ccPoly.m), zeros(obj.ccPoly.m,1), ...
                    obj.ccPoly.F', obj.K'*obj.sys.U.A(i,:)');
                L_HU(i,:) = lpSol';
            end
        end
        
        
        function P_RPI = getRPITerminalCost(obj)
            % we don't use 'full' option, so that P_RPI is symmetric (hence
            % we get positive semidefiniteness more easily)
            P_RPI = sdpvar(obj.ccPoly.m, obj.ccPoly.m);
            constr = [-obj.Q_RPI + P_RPI - obj.L_RPI'*P_RPI*obj.L_RPI >= 0;
                P_RPI >= 0];
            %             optimize(constr,trace(P_RPI),sdpsettings(obj.sdp_opts{:}));
            optimize(constr,trace(P_RPI),sdpsettings('solver','mosek','verbose',0));
            P_RPI=value(P_RPI);
        end
        
        function costMats = computeCostMats(obj, costMatricesSS)
            % from the constructor, we pass the cost matrices defined in
            % the state space (called here Q*_x/u,*=c,v). Due to comparison
            % reasons, we must redefine the system in the parameter space,
            % so transformations must be carried out. For this reason, we
            % save all the matrices in a single struct.
            
            % frequently used constant throughout function
            nx = obj.sys.nx; nu = obj.sys.nu;
            m = obj.ccPoly.m; m_bar = obj.ccPoly.m_bar;
            
            % define the cost struct
            costMats = struct;
            
            % CCPolytope volume cost (in state space)
            costMats.Qv_x = costMatricesSS.Qv(1:nx, 1:nx);
            costMats.Qv_u = costMatricesSS.Qv(nx+1:end, nx+1:end);
            
            % CCPolytope volume cost (in parameter space)
            V_bar = mean(cat(3, obj.ccPoly.Vi_s{:}), 3); % each Vi_s{i} is a 2D matrix, averaging along 3rd dim.
            U_bar = (1/m_bar)*repmat(eye(nu),1,m_bar); % nu x (nu*m_bar)
            
            Qv_yp = zeros(m,m);
            for i = 1:m_bar
                Qv_yp = Qv_yp + (obj.ccPoly.Vi_s{i}-V_bar)'*costMats.Qv_x*(obj.ccPoly.Vi_s{i}-V_bar);
            end
            costMats.Qv_yp = Qv_yp;
            
            Qv_up = zeros(nu*m_bar, nu*m_bar);
            Iu_mats = {};
            for i = 1:m_bar
                Iu_mats{i} = sparse(nu,nu*m_bar);
                Iu_mats{i}(:,(i-1)*nu+1:i*nu) = eye(nu);
                Qv_up = Qv_up + (Iu_mats{i}-U_bar)'*costMats.Qv_u*(Iu_mats{i}-U_bar);
            end
            costMats.Qv_up = Qv_up;
            
            % CCPolytope center cost (in state space)
            costMats.Qc_x = costMatricesSS.Qc(1:nx, 1:nx);
            costMats.Qc_u = costMatricesSS.Qc(nx+1:end, nx+1:end);
            
            % CCPolytope center cost (in parameter space)
            costMats.Qc_yp = V_bar'*costMats.Qc_x*V_bar;
            costMats.Qc_up = V_bar'*costMats.Qc_u*V_bar;
            
            % reference tracking (in state space)
            costMats.Qr = costMatricesSS.Qr;
        end
    end
    
end

