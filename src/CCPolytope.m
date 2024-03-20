classdef CCPolytope % TODO: does it need to be an handle class?
    %CCPOLYTOPE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        sys % dynamical system associated to the polytope
        ccPoly % Polyhedron
        m % number of facets
        m_bar % number of vertices
        Vi_s % cell array of vertices
        E % configuration-constrained polytope of parameter
        
        d % polytopic enclosure of the uncertainty set
        
    end
    
    properties(Dependent)
        F % in H-rep (Ax<=b) for the CCPolytope, F=A. Defined for brevity
    end
    
    methods % GETTER methods
        function F = get.F(obj)
            F = obj.ccPoly.A;
        end
    end
    
    methods
        function obj = CCPolytope(varargin)
            % Two admissible cases for inputs:
            % {dynSys, F} : use directly a provided CC-template (useful for dim=2)
            % {dynSys, C_tilde, D} : compute a valid CC-template (suggested for dim >=3)
            switch nargin
                case 2
                    F = varargin{2};
                case 3
                    F = obj.getCCTemplateMatrixCASADI(varargin{:});
                otherwise
                    error("Number of arguments provided : %d. Only 2 and 3 accepted.",length(varargin))
            end
            obj.sys = varargin{1};
            
            obj.m = size(F,1);
            y0 = ones(obj.m,1); % y0 = sigma
            
            obj.ccPoly = Polyhedron(F,y0);
            obj.Vi_s = obj.getVi_s(); % TODO: use a GETTER method
            obj.m_bar = length(obj.Vi_s);
            
            obj.E = sparse(obj.computeEMatrix());
            
            obj.d = obj.sys.W_dist.support(F');
        end
    end
    
    methods (Access = private)
        
        function E = computeEMatrix(obj)
            E = zeros(obj.m_bar*obj.m, obj.m);
            for i=1:obj.m_bar
                E((i-1)*obj.m+1:i*obj.m,:) = obj.ccPoly.A*obj.Vi_s{i}-eye(obj.m);
            end
            E(abs(E)<1e-10)=0; % round numerical errors ~=0
        end
        
        function V = getVi_s(obj) % no longer needed the other argument
            vertices = num2cell(obj.ccPoly.V',1);
            V = cell(size(vertices));
            
            for i=1:length(vertices)
                % check for each individual vertex which constraint is active
                Vi_mask = abs(obj.ccPoly.A*vertices{i}-obj.ccPoly.b) <= 1e-10;
                
                % compute 1_Vi to get Vi=inv(Y)*1_Vi
                one_mat = zeros(obj.sys.nx,obj.m); % nx x m
                one_mat(sub2ind(size(one_mat),1:obj.sys.nx,find(Vi_mask)')) = 1;
                
                V{i} = obj.ccPoly.A(Vi_mask,:) \ one_mat;
            end
            
        end
        
        function F = getCCTemplateMatrixCASADI(~,sys,C_tilde,D)
            % - firstly, we get a user defined C_tilde (by using m and nx, we
            % decide template complexity)
            % - then, we compute with the opt. problem a proper transformation matrix
            % W_inv such that F = C_tilde*W_inv.
            nx = sys.nx;
            nu = sys.nu;
            
            Z_set = Polyhedron(C_tilde,ones(size(C_tilde,1),1) );
            zj_s =  num2cell(Z_set.V',1);
            
            %  optimization problem
            opti = casadi.Opti();
            
            W = opti.variable(nx,nx);
            M = opti.variable(nx,nx);
            uj_s = opti.variable(nu,length(zj_s));
            
            % distance definition
            X_V = num2cell(sys.X.V',1);
            
            epsilon =  opti.variable(size(D,1),1);
            
            d =  opti.variable(nx, length(X_V));
            s =  opti.variable(nx, length(X_V));
            sum_var = opti.variable(nx, length(X_V));
            
            for i=1:length(X_V) % all the points must fulfill the inclusion
                % minkowski sum
                opti.subject_to( sum_var(:,i) == d(:,i) + s(:,i) );
                opti.subject_to( D*d(:,i) <= epsilon );
                opti.subject_to( C_tilde*M*s(:,i) <= ones(size(C_tilde,1),1) );
                
                % inclusion
                opti.subject_to( X_V{i} == sum_var(:,i) );
            end
            
            % vertex-control problem
            % TODO: use lagrangian for disturbance.
            wl_s =  num2cell(sys.W_dist.V',1);
                        
            for j=1:length(zj_s) % vertices
                for i=1:length(sys.A_convh) % #models
                    for l=1:length(wl_s) % disturbance vertices
                        opti.subject_to( ...
                            C_tilde*M*(sys.A_convh{i}*W*zj_s{j}+sys.B_convh{i}*uj_s(:,j)+wl_s{l})<=ones(size(C_tilde,1),1) ...
                            );
                    end
                end
                opti.subject_to( sys.X.A*W*zj_s{j} <= sys.X.b );
                opti.subject_to( sys.U.A*uj_s(:,j) <= sys.U.b );
            end
            opti.subject_to( W*M == eye(nx) );
            
            costFun = sum(abs(epsilon));
            
            opti.minimize(costFun);
            
            % set print_level=0 and print_time=0 to remove any debug output
            ipopt = struct('print_level',3,'sb','yes');
            opti_opts = struct('ipopt', ipopt, 'print_time', 1);
            opti.solver('ipopt',opti_opts)
            
            % non-zero initial condition for matrices
            opti.set_initial(W, eye(nx));
            opti.set_initial(M,	eye(nx));
            
            sol = opti.solve();
            
            F = C_tilde*sol.value(M);
        end
    end    
end

