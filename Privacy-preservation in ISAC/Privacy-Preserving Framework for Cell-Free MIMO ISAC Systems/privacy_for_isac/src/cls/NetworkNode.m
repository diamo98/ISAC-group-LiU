classdef NetworkNode < CommSystemObj
    %NETWORKNODE Summary of this class goes here
    %   Detailed explanation goes here'

    properties (GetAccess=public, SetAccess=protected)
        nAntennas(1, 1) {mustBeNaturalNumber(nAntennas)} = 0 % Number of antennas
    end
  
    properties (Access=public)
        precoderMatrix = 0;
        decoderMatrix = 0;
    end
    
    methods
        function obj = NetworkNode(pos, nAntennas)
            %NETWORKNODE Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@CommSystemObj(pos);
            obj.nAntennas = nAntennas;
        end
        
        % TODO: 
        % dir?
        function a = calcSteeringVector(obj, argin, dir)
            %CALCSTEERINGVECTOR
            %   Assumed half wavelength equally spaced linear antenna array
            if isa(argin, "CommSystemObj")
                if strcmp(dir, "from")
                    theta = obj.angleTo(argin) + pi;
                elseif strcmp(dir, "to")
                    theta = obj.angleTo(argin);
                else
                    error("not valid dir, must be from or to");
                end
                a = [exp(1j*pi*(0:(obj.nAntennas-1))*cos(theta))].';
            elseif isa(argin, "double")
                a = [exp(1j*pi*(0:(obj.nAntennas-1))*cos(argin))].';
            else
                disp(argin)
                error("Bad datatype in calcSterringVector");
            end
        end
    end
end

