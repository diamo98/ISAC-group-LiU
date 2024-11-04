classdef (Abstract) CommSystemObj < handle
    %COMMSYSTEMOBJ Base object class for the modelled communication system

    properties (GetAccess=public, SetAccess=private)
        pos(2, 1) {validateattributes(pos, {'numeric'}, {'finite', 'real'})} % Position of object
    end
    
    methods
        function obj = CommSystemObj(pos)
            %COMMSYSTEMOBJ
            % Arguments:
            % pos: Position of object
            p = inputParser;
            addRequired(p, 'pos');
            parse(p, pos);
            obj.pos = reshape(p.Results.pos, 2, []);
        end
    end

    methods (Sealed, Access=public)
        function dist = distanceTo(obj, other)
            %DISTANCETO
            % Arguments:
            % other: CommSystemObj to calculate distance to
            arguments
                obj
                other(1,1) CommSystemObj
            end

            dist = norm(other.pos-obj.pos);
        end

        function dir = angleTo(obj, other)
            %ANGLETO
            % Arguments:
            % other: CommSystemObj to calculate angle to
            arguments
                obj
                other(1,1) CommSystemObj
            end
            dir = angle(other.pos(1)-obj.pos(1) + (other.pos(2)-obj.pos(2))*1j);
        end
    end
end
