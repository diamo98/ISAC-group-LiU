classdef AccessPoint < NetworkNode
    %ACCESSPOINT Summary of this class goes here
    %   Detailed explanation goes here
    methods
        function obj = AccessPoint(pos, nAntennas)
            %ACCESSPOINT Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@NetworkNode(pos, nAntennas);
        end
    end
end
