classdef DynSystem < handle
    %DYNSYS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        A_convh %
        B_convh
        X
        U
        W_dist
    end
    
    properties(Dependent)
        nx
        nu
    end
    
    methods 
        function nx = get.nx(obj)
            nx = size(obj.A_convh{1},2);
        end
        function nu = get.nu(obj)
            nu = size(obj.B_convh{1},2);
        end
    end
    
    methods (Access = public)
        function obj = DynSystem(A_convh, B_convh, X, U, W_dist)
            %DYNSYS Construct an instance of this class
            %   Detailed explanation goes here
            obj.A_convh = A_convh;
            obj.B_convh = B_convh;
            obj.X = X;
            obj.U = U;
            obj.W_dist = W_dist;
        end
        
        function updateSysMatrices(obj, A_convh_new, B_convh_new)
        % used for online model refinement
            obj.A_convh = A_convh_new;
            obj.B_convh = B_convh_new;
        end
    end
end

