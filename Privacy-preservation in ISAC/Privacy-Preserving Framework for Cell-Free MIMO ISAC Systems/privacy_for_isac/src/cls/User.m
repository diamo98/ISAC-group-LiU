classdef User < NetworkNode
    %USER Summary of this class goes here
    %   Detailed explanation goes here
    methods
        function obj = User(pos, nAntennas)
            %USER Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@NetworkNode(pos, nAntennas);
        end
    end
end

