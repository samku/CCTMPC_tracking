classdef CCPolytope
    
    properties (SetAccess = private)
        sys % dynamical system associated to the polytope
        ccpolytope % Configuration Constrained polytope
        m % number of facets
        m_bar % number of vertices
        Vi_s % cell array of vertices
        E % configuration-constrained polytope of parameter
        
        d % support (polytopic enclosure) of the uncertainty set
    end
    
    properties(Dependent)
        F % in H-rep (Ax<=b) for the CCPolytope, F=A. Defined for brevity
    end
    
    methods % GETTER methods
        function F = get.F(obj)
            F = obj.ccpolytope.A;
        end
    end
    
    methods (Access = public)
        function obj = CCPolytope(sys, varargin)
            % Two admissible constructor signatures:
            % {sys, F} : use a provided CC-template directly (useful for dim=2)
            % {sys, C_tilde, D} : compute a valid CC-template (suggested for dim >=3)
            obj.sys = sys;
            
            switch length(varargin)
                case 1
                    F = varargin{1};
                case 2
                    F = obj.getCCTemplateMatrixCASADI(varargin{:});
                otherwise
                    error("Number of arguments provided : %d. Only 2 and 3 accepted.",nargin)
            end
            obj.m = size(F,1);
            
            % using ones() as initial parameter, the CCPolytope is entirely simple.
            obj.ccpolytope = Polyhedron(F,ones(obj.m,1));
            
            obj.Vi_s = obj.getVi_s();
            obj.m_bar = length(obj.Vi_s);
            
            obj.E = sparse(obj.computeEMatrix()); % usually very sparse.
            
            obj.d = obj.sys.W_dist.support(F');
        end
    end
    
    methods (Access = private)
        
        function E = computeEMatrix(obj)
            % conic constraint defining the vertex configuration domain explicitly
            E = zeros(obj.m_bar*obj.m, obj.m);
            for i=1:obj.m_bar
                E((i-1)*obj.m+1:i*obj.m,:) = obj.F*obj.Vi_s{i}-eye(obj.m);
            end
            E(abs(E)<1e-10) = 0; % round numerical errors ~=0
        end
        
        
        function V = getVi_s(obj)
            % computing Vi_s such that Vi*y0 is an isolated vertex of the CCPolytope
            vertices = num2cell(obj.ccpolytope.V',1);
            V = cell(size(vertices));
            
            for i=1:length(vertices)
                % check for each individual vertex which constraint is active
                Vi_mask = abs(obj.F*vertices{i}-obj.ccpolytope.b) <= 1e-10;
                
                % compute 1_Vi to get Vi=inv(Y)*1_Vi
                one_mat = zeros(obj.sys.nx,obj.m); % nx x m
                one_mat(sub2ind(size(one_mat),1:obj.sys.nx,find(Vi_mask)')) = 1;
                
                V{i} = obj.F(Vi_mask,:) \ one_mat;
            end
        end
        
        
        function F = getCCTemplateMatrixCASADI(obj,C_tilde,D)
            % from [https://arxiv.org/abs/2309.02384,Section IV-D]:
            % - firstly, we define a C_tilde matrix (m,nx), thus fixing
            % the polytope complexity)
            % - then, we solve an opt. problem to find the proper
            % transformation matrix W_inv such that F = C_tilde*W_inv is a
            % template for a CCPolytope.
            nx = obj.sys.nx; nu = obj.sys.nu;
            
            Z_set = Polyhedron(C_tilde,ones(size(C_tilde,1),1) );
            zj_s = num2cell(Z_set.V',1);
            
            %  optimization problem
            opti = casadi.Opti();
            
            W = opti.variable(nx,nx);
            M = opti.variable(nx,nx);
            uj_s = opti.variable(nu,length(zj_s));
            
            % distance definition
            X_V = num2cell(obj.sys.X.V',1);
            
            epsilon =  opti.variable(size(D,1),1);
            
            d_var =  opti.variable(nx, length(X_V));
            s_var =  opti.variable(nx, length(X_V));
            sum_var = opti.variable(nx, length(X_V));
            
            for i=1:length(X_V) % all the points must fulfill the inclusion
                % minkowski sum
                opti.subject_to( sum_var(:,i) == d_var(:,i) + s_var(:,i) );
                opti.subject_to( D*d_var(:,i) <= epsilon );
                opti.subject_to( C_tilde*M*s_var(:,i) <= ones(size(C_tilde,1),1) );
                
                % inclusion
                opti.subject_to( X_V{i} == sum_var(:,i) );
            end
            
            % vertex-control problem
            % TODO: use multipliers instead of vertices for disturbance
            wl_s =  num2cell(obj.sys.W_dist.V',1);
            
            for j=1:length(zj_s) % vertices
                for i=1:length(obj.sys.A_convh) % #models
                    for l=1:length(wl_s) % disturbance vertices
                        opti.subject_to( ...
                            C_tilde*M*(obj.sys.A_convh{i}*W*zj_s{j}+obj.sys.B_convh{i}*uj_s(:,j)+wl_s{l})<=ones(size(C_tilde,1),1) ...
                            );
                    end
                end
                opti.subject_to( obj.sys.X.A*W*zj_s{j} <= obj.sys.X.b );
                opti.subject_to( obj.sys.U.A*uj_s(:,j) <= obj.sys.U.b );
            end
            opti.subject_to( W*M == eye(nx) );
            
            costFun = sum(abs(epsilon)); % 1-norm
            
            opti.minimize(costFun);
            
            % set print_level=0 and print_time=0 to remove any debug output
            ipopt = struct('print_level',3,'sb','yes');
            opti_opts = struct('ipopt', ipopt, 'print_time', 1);
            opti.solver('ipopt',opti_opts)
            
            % non-zero initial condition for matrices (improving OCP performance)
            opti.set_initial(W, eye(nx));
            opti.set_initial(M,	eye(nx));
            
            sol = opti.solve();
            
            F = C_tilde*sol.value(M);
        end
        
    end
    
end