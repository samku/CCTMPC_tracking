classdef CostFunctionManager
    % (?) HANDLE not required, just a data class.
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        cost_RCI % must be convex function (not necessarily strictly)
        cost_L_k % strictly convex function (Q > 0 )
        cost_L_N % strictly convex function (P > 0 && Q <= (1-gamma^2)*P)
    end
    
    methods
        function obj = CostFunctionManager(sys, ccpoly, RCI_cost, Stage_cost, Term_cost)
            % if the value is a struct of matrices, then construct the cost
            % functions; otherwise, set the function handle accordingly
            
            switch class(RCI_cost)
                case 'struct'
                    obj.cost_RCI = obj.buildRCICost(sys, ccpoly, RCI_cost);
                case 'function_handle'
                    obj.cost_RCI = RCI_cost;
                otherwise
                    error("Unexpected DataType: %s. Only 'struct' and 'function_handle' accepted.",class(RCI_cost))
            end
            
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            switch class(Stage_cost)
                case 'double'
                    % TODO: should we assert (Q > 0)?
                    % TODO: assert dimensions! 
                    obj.cost_L_k = @(y_k,u_k,ys,us) weightSq2Norm([y_k-ys; u_k(:)-us(:)],Stage_cost);
                case 'function_handle'
                    obj.cost_L_k = Stage_cost;
                otherwise
                    error("Unexpected DataType: %s. Only 'struct' and 'function_handle' accepted.",class(Stage_cost))
            end
            
            switch class(Term_cost)
                case 'double'
                    % TODO: should we assert (P > 0 && Q <= (1-gamma^2)*P)?
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
            V_bar = mean(cat(3, ccpoly.Vi_s{:}), 3); % each Vi_s{i} is a 2D matrix, concatenating on 3rd dim.
            U_bar = (1/ccpoly.m_bar)*repmat(eye(sys.nu),1,ccpoly.m_bar); % nu x (nu*m_bar)
            
            A_bar = mean(cat(3, sys.A_convh{:}), 3);
            B_bar = mean(cat(3, sys.B_convh{:}), 3);
            
            M = null([A_bar-eye(sys.nx) B_bar]);
            
            weightSq2Norm = @(vector, mat) vector'*mat*vector;
            
            cost_RCI = @(ys, us, r, th_s)...        
                arrayfun(@(i)sum( weightSq2Norm(  [(V_bar-ccpoly.Vi_s{i})*ys ; U_bar*us(:)-us(:,i)], RCI_cost.Qv)  ), 1:ccpoly.m_bar ) + ...
                weightSq2Norm( blkdiag(V_bar,U_bar)*[ys;us(:)] - M*th_s, RCI_cost.Qc) + ...
                weightSq2Norm( r - [sys.C, sys.D]*M*th_s, RCI_cost.Qr);

        end        
    end
end













