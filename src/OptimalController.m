classdef OptimalController < handle
    %OPTIMALCONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dynSys
        N
    end
    
    methods
        function obj = OptimalController(dynSys, N)
            %OPTIMALCONTROLLER Construct an instance of this class
            %   Detailed explanation goes here
            obj.dynSys = dynSys;
            obj.N = N;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end