classdef CostFunctionManager
    % Collector class to handle objective function definition, that must be
    % later passed to the OptimalController class.
    
    properties (SetAccess = private)
        cost_RCI % must be a convex function (not necessarily strictly)
        cost_L_k % strictly convex function (for quadratic cost, Q > 0)
        cost_L_N % strictly convex function (for quadratic cost, P > 0 && Q <= (1-gamma^2)*P)
    end
    
    methods (Access = public)
        function obj = CostFunctionManager(sys, ccPoly, RCI_cost, Stage_cost, Term_cost)
            % if matrices are passed, construct a quadratic cost function;
            % otherwise, set the function handle accordingly.
            % NOTE: in the second case, function signature must be:
            % RCI_cost:     @(ys,us,r,th_s)
            % Stage_cost:   @(y_k,u_k,ys,us)
            % Term_cost:    @(y_N,u_N,ys,us)
            
            switch class(RCI_cost)
                case 'struct'
                    obj.cost_RCI = obj.buildRCICost(sys, ccPoly, RCI_cost);
                case 'function_handle'
                    obj.cost_RCI = RCI_cost;
                otherwise
                    error("Unexpected DataType: %s. Only 'struct' and 'function_handle' accepted.",class(RCI_cost))
            end
            
            weightSq2Norm = @(vector, mat) vector'*mat*vector; % squared 2-norm
            switch class(Stage_cost)
                case 'double'
                    obj.cost_L_k = @(y_k,u_k,ys,us) weightSq2Norm([y_k-ys; u_k(:)-us(:)],Stage_cost);
                case 'function_handle'
                    obj.cost_L_k = Stage_cost;
                otherwise
                    error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(Stage_cost))
            end
            
            switch class(Term_cost)
                case 'double'
                    obj.cost_L_N = @(y_N,u_N,ys,us) weightSq2Norm([y_N-ys; u_N(:)-us(:)],Term_cost);
                case 'function_handle'
                    obj.cost_L_N = Term_cost;
                otherwise
                    error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(RCI_cost))
            end
        end
    end
    
    methods (Access = private)
        
        function cost_RCI = buildRCICost(~,sys,ccpoly,RCI_cost)
            % refer to Section II-D from the paper.
            
            V_bar = mean(cat(3, ccpoly.Vi_s{:}), 3); % each Vi_s{i} is a 2D matrix, averaging along 3rd dim.
            U_bar = (1/ccpoly.v)*repmat(eye(sys.nu),1,ccpoly.v); % nu x (nu*v)
            
            A_bar = mean(cat(3, sys.A_convh{:}), 3);
            B_bar = mean(cat(3, sys.B_convh{:}), 3);
            
            M = null([A_bar-eye(sys.nx) B_bar]);
            
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            
            cost_RCI = @(ys, us, r, th_s)...
                arrayfun(@(i)sum( weightSq2Norm(  [(V_bar-ccpoly.Vi_s{i})*ys ; U_bar*us(:)-us(:,i)], RCI_cost.Qv)  ), 1:ccpoly.v ) + ...
                weightSq2Norm( blkdiag(V_bar,U_bar)*[ys;us(:)] - M*th_s, RCI_cost.Qc) + ...
                weightSq2Norm( r - [sys.C, sys.D]*M*th_s, RCI_cost.Qr);
        end
    end
end













