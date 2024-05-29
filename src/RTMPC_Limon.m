classdef RTMPC_Limon < handle
    % Rigid Tube MPC formulation for set-point tracking, based on Limon's 
    % paper implementation. For this reason, it can be used only with LTI 
    % systems. Differently from the original formulation, we also use the
    % parameter space (and not only the state space), giving us a fairer
    % comparison with our method.
    
    properties (SetAccess = private)        
        optContrSol % struct, solution from OCP
        
        feasRegion % Polytope, (approximated) feasible region for the OCP
        
        sys % dynamical system 
        ccPoly % vertex and parameters from the CCPolytope are used for comparison
        costMats % struct that aggregates all the relevant matrices for the system
        N % prediction horizon for the OCP
        lambda % computing a tight approx. of the mRPI set as terminal set
        sdp_opts % YALMIP optimizer options (solver, verbosity...)
        
        y_mRPI, u_mRPI
        termSetAug, hX_tight, hU_tight
        M, P
    end
    
    properties (Access = private)
        optController % YALMIP optimizer, Optimal Controller
        mpcLaw % YALMIP optimizer, compute the input in state (and not parameter) space
        
        directionalFeasRegion % YALMIP optimizer, compute the farthest feasible
        % initial state with respect to a specific facet.
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
        function obj = RTMPC_Limon(sys, ccPoly, costMatricesSS, N, lambda, sdp_opts)
            assert(sys.np == 1, "Cannot proceed. RTMPC from Limon is defined for only LTI systems.")
            
            obj.sys = sys;
            obj.ccPoly = ccPoly;
            obj.N = N;
            obj.lambda = lambda;
            obj.sdp_opts = sdp_opts;
            
            % struct{Qv_(y/u)p, Qv_(x/u), Qc_(y/u)p, Qc_(x/u), Qr}
            obj.costMats = obj.computeCostMats(costMatricesSS); %{Qv, Qc, Qr}
            
            % compute terminal set for [z;theta]. Note that some class
            % attributes are instantiated inside the method.
            obj.termSetAug = obj.computeTerminalSet();
            
            % instantiate YALMIP optimizers
            obj.optController = obj.initOptController();
            obj.mpcLaw = obj.initMPCLaw();
            
            obj.optContrSol = cell(4,1); % {y,u,theta,cost} fields

            obj.directionalFeasRegion = obj.initDirFeasRegion();
            obj.feasRegion = nan; % computation (slow) delayed, as only used in plots
        end
        
        
        function uMPC = solve(obj, x_curr, r_curr)
            % solve OCP (solution in parameter space)
            obj.optContrSol = obj.optController(x_curr, r_curr);
            uMPC = obj.mpcLaw(x_curr, obj.optContrSol{1}(:,1), obj.optContrSol{2}(:,1));
            
            if isnan(uMPC)
                if sum(abs(obj.findFeasibleX0(x_curr,eye(obj.sys.nx)) - x_curr)) > 1e-5
                    error("The current state does not lie inside the feasible"+...
                        " region of the control scheme. Consider changing initial condition.")
                end
                error("uMPC returned NaN. Consider enabling debug mode.")
            end
        end
        
        
        function feas_x0 = findFeasibleX0(obj, x0_des, costFun)
            % find the closest feasible initial state to the desired one
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            Mx = obj.M(1:obj.sys.nx,:); Mu = obj.M(obj.sys.nx+1:end,:);
            
            % optContr YALMIP definitions
            theta = sdpvar(obj.sys.nu, 1,'full');
            x0 = sdpvar(obj.sys.nx, 1, 'full');
            x = sdpvar(obj.sys.nx, obj.N,'full');
            u = sdpvar(obj.sys.nu, obj.N-1,'full');
            
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
            constr = [obj.ccPoly.F*x0 <= obj.y_mRPI + obj.ccPoly.F*x(:,1); % initial condition
                obj.sys.X.A*Mx*theta <= obj.sys.X.b;
                obj.sys.U.A*Mu*theta <= obj.sys.U.b];
            
            % stage cost and constraints
            for t = 1:obj.N-1
                constr=[constr; x(:,t+1) == obj.sys.A_curr()*x(:,t)+obj.sys.B_curr()*u(:,t)];
                constr=[constr; obj.sys.X.A*x(:,t) <= obj.hX_tight];
                constr=[constr; obj.sys.U.A*u(:,t) <= obj.hU_tight];
            end
            % terminal constraints
            constr = [constr; obj.termSetAug.A*[x(:,obj.N); theta] <= obj.termSetAug.b];
            
            % solve the OCP
            optimize(constr,cost,sdpsettings(obj.sdp_opts{:}));
            feas_x0 = value(x0);
        end
        
        
        function feasRegion = computeFeasRegion(obj, n_facets)
            % compute a new feasible region given a template of #facets complexity.
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
            
            % compute mRPI facets from CCPolytope
            [obj.y_mRPI, obj.u_mRPI] = obj.get_mRPI_facets();

            % compute the steady-state subspace
            obj.M = null([A-eye(nx) B]);
            Mx = obj.M(1:nx,:); Mu = obj.M(nx+1:end,:);
            
            % get affine feedback law from LQR solution
            [K, obj.P] = dlqr(A, B, obj.costMats.Qc_x, obj.costMats.Qc_u);
            K = -K;
            
            % extend state vector to [z;theta]
            A_aug = [A+B*K, B*Mu-B*K*Mx; 
                zeros(nu,nx), eye(nu)];
            
            % compute tightened X/U constraint sets
            obj.hX_tight = zeros(size(hX)); obj.hU_tight = zeros(size(hU));
            for i=1:length(hX)
                xx = cplexlp(-HX(i,:)',obj.ccPoly.F,obj.y_mRPI);
                obj.hX_tight(i) = hX(i) - HX(i,:)*xx;
            end
            U_mRPI_set = Polyhedron(obj.u_mRPI').minHRep();
            HU_mRPI = U_mRPI_set.A; hU_mRPI = U_mRPI_set.b;
            for i=1:length(hU)
                uu = cplexlp(-HU(i,:)',HU_mRPI,hU_mRPI);
                obj.hU_tight(i) = hU(i) - HU(i,:)*uu;
            end
            
            % Compute the augmented constraint set
            H_aug = [HX, zeros(size(HX,1),nu);
                HU*K, HU*Mu-HU*K*Mx;
                zeros(size(HX,1),nx), HX*Mx;
                zeros(size(HU,1),nx), HU*Mu];
            
            h_aug = [obj.hX_tight;
                obj.hU_tight;
                obj.lambda*obj.hX_tight;
                obj.lambda*obj.hU_tight];
            
            % iterate constraints through system dynamics to overapproximate
            % the terminal set from above (volume of the set diminishes at
            % each iteration)
            H_RPI = H_aug; h_RPI = h_aug;
            for t=1:100 % hard-coded exit condition, usually suff. tight approx.
                H_LHS = [H_RPI*A_aug; H_aug];   
                h_LHS = [h_RPI; h_aug];
                termSetAug = Polyhedron(H_LHS,h_LHS).minHRep();
                H_RPI = termSetAug.A;    
                h_RPI = termSetAug.b;
                
                % Verify invariance H(Ax+d)<=h for all x:Hx<=h
                new_h_RPI = zeros(size(h_RPI));
                for i=1:size(h_RPI,1)
                    xx = cplexlp(-(H_RPI(i,:)*A_aug)',H_RPI,h_RPI);
                    new_h_RPI(i) = H_RPI(i,:)*A_aug*xx;
                end
                if max(new_h_RPI-h_RPI)<=1e-6
                   break
                end
            end
        end
        
        
        function optContr = initOptController(obj)            
            % frequently used constants throughout function
            A = obj.sys.A_curr(); B = obj.sys.B_curr(); C = obj.sys.C; D = obj.sys.D;
            
            Mx = obj.M(1:obj.sys.nx,:); Mu = obj.M(obj.sys.nx+1:end,:);
                        
            % solve OCP
            r = sdpvar(obj.sys.ny, 1,'full');
            theta = sdpvar(obj.sys.nu, 1,'full');
            x0 = sdpvar(obj.sys.nx, 1, 'full');
            x = sdpvar(obj.sys.nx, obj.N,'full');
            u = sdpvar(obj.sys.nu, obj.N-1,'full');
            
            % initial condition & steady-state constraints
            constr = [obj.ccPoly.F*x0 <= obj.y_mRPI + obj.ccPoly.F*x(:,1); % initial condition
                obj.sys.X.A*Mx*theta <= obj.sys.X.b;
                obj.sys.U.A*Mu*theta <= obj.sys.U.b];
            % tracking cost
            cost = (r - [C, D]*obj.M*theta)'*obj.costMats.Qr*(r - [C, D]*obj.M*theta);
            
            % stage cost and constraints
            for t = 1:obj.N-1
                constr=[constr; x(:,t+1) == A*x(:,t)+B*u(:,t)];
                constr=[constr; obj.sys.X.A*x(:,t) <= obj.hX_tight];
                constr=[constr; obj.sys.U.A*u(:,t) <= obj.hU_tight];
                
                cost = cost + (x(:,t)-Mx*theta)'*obj.costMats.Qc_x*(x(:,t)-Mx*theta)...
                    + (u(:,t)-Mu*theta)'*obj.costMats.Qc_u*(u(:,t)-Mu*theta);
            end
            % terminal cost and constraints
            constr = [constr; obj.termSetAug.A*[x(:,obj.N); theta] <= obj.termSetAug.b];
            cost = cost + (x(:,obj.N)-Mx*theta)'*obj.P*(x(:,obj.N)-Mx*theta);
            
            % instantiate Optimizer
            optContr = optimizer(constr,cost,sdpsettings(obj.sdp_opts{:}),...
                {x0,r},...      % input parameters
                {x,u,theta,cost});  % output parameters
        end
        
        
        function mpcLaw = initMPCLaw(obj)
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            
            x0Opt = sdpvar(obj.sys.nx, 1,'full'); % centroid of CCpolytope (OCP solution)
            u0Opt = sdpvar(obj.sys.nu, 1,'full'); % open-loop control component (OCP solution)
            x_t = sdpvar(obj.sys.nx, 1,'full'); % current state of the system
            
            thetas = sdpvar(obj.ccPoly.v, 1,'full');
            
            constr = [];
            constr = [constr; sum(thetas)==1; thetas >= 0];
            
            x_vert = [];
            for k=1:obj.ccPoly.v
                x_vert = [x_vert, obj.ccPoly.Vi_s{k}*(obj.y_mRPI+obj.ccPoly.F*x0Opt)];
            end
            
            u_vert = repmat(u0Opt,1,obj.ccPoly.v) + obj.u_mRPI; 
            
            constr = [constr; x_t == x_vert*thetas];
            
            
            cost = weightSq2Norm(thetas, eye(size(thetas,1))); % squared 2-norm
            
            mpcLaw = optimizer(constr, cost, sdpsettings(obj.sdp_opts{:}), ...
                {x_t, x0Opt, u0Opt}, ... % input parameters
                {u_vert*thetas} ... % output parameters
                ); % (nu x v) * (v x 1) -> (nu x 1)
        end
        
        
        function dirFeasRegion = initDirFeasRegion(obj)
            % compute the directional feasibility region with respect to a
            % specific direction.
            Mx = obj.M(1:obj.sys.nx,:); Mu = obj.M(obj.sys.nx+1:end,:);
            c_vec = sdpvar(obj.sys.nx,1,'full');
            
            theta = sdpvar(obj.sys.nu,1,'full');
            x0 = sdpvar(obj.sys.nx,1,'full');
            x = sdpvar(obj.sys.nx,obj.N,'full');
            u = sdpvar(obj.sys.nu,obj.N-1,'full');
            
            % cost function
            cost = -c_vec'*x0; % maximize the direction
            
            % initial condition & steady-state constraints
            constr = [obj.ccPoly.F*(x0-x(:,1)) <= obj.y_mRPI; % initial condition
                obj.sys.X.A*Mx*theta <= obj.sys.X.b;
                obj.sys.U.A*Mu*theta <= obj.sys.U.b];
            
            % stage constraints
            for t = 1:obj.N-1
                constr=[constr; x(:,t+1) == obj.sys.A_curr()*x(:,t)+obj.sys.B_curr()*u(:,t)];
                constr=[constr; obj.sys.X.A*x(:,t) <= obj.hX_tight];
                constr=[constr; obj.sys.U.A*u(:,t) <= obj.hU_tight];
            end
            % terminal constraints
            constr = [constr; obj.termSetAug.A*[x(:,obj.N); theta] <= obj.termSetAug.b];
            
            dirFeasRegion = optimizer(constr, cost, sdpsettings(obj.sdp_opts{:}),...
                {c_vec},... % input parameters
                {x0});      % output parameters
        end
        
        
        function [y_mRPI, u_mRPI] = get_mRPI_facets(obj)
            % Limon mRPI set
            A = obj.sys.A_curr(); B = obj.sys.B_curr();
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            
            ys = sdpvar(obj.ccPoly.m, 1,'full');
            us = sdpvar(obj.sys.nu, obj.ccPoly.v,'full');
            
            % define constraints
            constr = [];
            constr = [constr; obj.ccPoly.E*ys <= 0];
            
            for k = 1:obj.ccPoly.v
                % Robust Forward Invariant constraint
                constr = [constr; obj.ccPoly.F*(A*obj.ccPoly.Vi_s{k}*ys + B*us(:,k)) + obj.ccPoly.d <= ys];
                % state and input constraints
                constr = [constr; obj.sys.X.A*obj.ccPoly.Vi_s{k}*ys <= obj.sys.X.b;   obj.sys.U.A*us(:,k)<= obj.sys.U.b];
            end
            
            % define cost function
            cost = weightSq2Norm(ys, obj.costMats.Qv_yp+obj.costMats.Qc_yp) ...
                + weightSq2Norm(us(:), obj.costMats.Qv_up+obj.costMats.Qc_up);
            
            % solve the OCP
            optimize(constr,cost,sdpsettings(obj.sdp_opts{:}));
            
            y_mRPI = value(ys);
            u_mRPI = value(us);
        end
        
        
        function costMats = computeCostMats(obj, costMatricesSS)
            % from the constructor, we pass the cost matrices defined in
            % the state space (called here Q*_x/u,*=c,v). For comparison
            % purposes, we redefine the system in the parameter space. All 
            % the matrices are saved in a single struct.
            
            % frequently used constant throughout function
            nx = obj.sys.nx; nu = obj.sys.nu;
            m = obj.ccPoly.m; v = obj.ccPoly.v;
            
            % define the cost struct
            costMats = struct;
            
            % CCPolytope volume cost (in state space)
            costMats.Qv_x = costMatricesSS.Qv(1:nx, 1:nx);
            costMats.Qv_u = costMatricesSS.Qv(nx+1:end, nx+1:end);
            
            % CCPolytope volume cost (in parameter space)
            V_bar = mean(cat(3, obj.ccPoly.Vi_s{:}), 3); % each Vi_s{i} is a 2D matrix, averaging along 3rd dim.
            U_bar = (1/v)*repmat(eye(nu),1,v); % nu x (nu*v)
            
            Qv_yp = zeros(m,m);
            for i = 1:v
                Qv_yp = Qv_yp + (obj.ccPoly.Vi_s{i}-V_bar)'*costMats.Qv_x*(obj.ccPoly.Vi_s{i}-V_bar);
            end
            costMats.Qv_yp = Qv_yp;
            
            Qv_up = zeros(nu*v, nu*v);
            Iu_mats = {};
            for i = 1:v
                Iu_mats{i} = sparse(nu,nu*v);
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