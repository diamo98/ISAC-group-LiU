classdef Target < CommSystemObj
    %TARGET Summary of this class goes here
    %   Detailed explanation goes here

    % Property and argument validation!
    properties (SetAccess=private, GetAccess=public)
        rcs = 0;
    end

    methods
        function obj = Target(pos, rcs)
            %TARGET Construct an instance of this class
            %   Detailed explanation goes here
            %obj.Property1 = inputArg1 + inputArg2;
            obj = obj@CommSystemObj(pos)
            obj.rcs = rcs;
        end
    end
end

