classdef SimulationEnvironment < handle & matlab.mixin.Copyable
    %SIMULATIONENVIRONMENT Summary of this class goes here
    %   Detailed explanation goes here

    % TODO: Gå igenom metoderna så alla är kallningsbara i vilken ordning
    % som helst (med varierande värdefullt resultat förstås)
    
    % SET ACCESS = PROTECTED
    properties (SetAccess=public, GetAccess=public)
        % TODO:
        % Property validation!

        nObservations = 0;
        noiseVar = 0;

        % Communication
        cSignal = 0;
        cH = 0;
        cBetas = 0;

        % Sensing
        sSignal = 0;
        sH = 0;
        sBetas = 0;

        % Communication System Objects
        transmitAps = AccessPoint([0,0], 0);
        receiverAps = AccessPoint([0,0], 0);

        users = User([0,0], 0);
        adversaryIdx = 0;

        target = Target([0,0], 0);

        % Precoders
        %TODO: UNUSED??????????
        precoders = 0;

        % Covariance estimated
        estimatedR = cell(1,1);

        correctAngles = 0;
        estimatedAngles = 0;
    end
    
    methods
        function obj = SimulationEnvironment(nAntennas, userAntennas, ...
                apsLocs, receiverIdxs, usersLocs, adversaryIdx, targetLoc, ...
                rcsVar, noiseVar, maxPower, nObservations, seed)
            % input validation, like no more receivers than aps obviously and that there is a target
            % and locations and users and stuff!

            %SIMULATIONENVIRONMENT Construct an instance of this class
            %   Detailed explanation goes here
            % OBS: Create a copy of a SimulationEnvironment by calling the
            % constructor with another SimulationEnvironment as its
            % argument
            
            %TODO: bad praxis to not use varargin with names instead of using
            % nAntennas as a general input argument 
            if nargin == 1 && isa(nAntennas, 'SimulationEnvironment')
                obj = copy(nAntennas);
                return
            end

            % temp lösning
            if exist("seed", "var")
                rng(seed);
            end

            obj.noiseVar = noiseVar;
            obj.nObservations = nObservations;
            obj.adversaryIdx = adversaryIdx;

            for i = 1:height(usersLocs)
                obj.users(i) = User(usersLocs(i,:), userAntennas);
            end

            nUsers = numel(obj.users);

            for i = 1:height(apsLocs)
                % variance?
                accessPoints(i) = AccessPoint(apsLocs(i,:), nAntennas);

                % TODO:
                % default precoder is equal power and phase on all antenna
                % elements (don't use except for precoder comparisons)
                % Calculate the dBm instead of hardcode 100
                accessPoints(i).precoderMatrix = sqrt(maxPower/(nUsers+1))*ones(nAntennas, nUsers+1);
            end
            apIdxs = logical(zeros(1, numel(accessPoints)));
            apIdxs(receiverIdxs) = 1;
            obj.transmitAps = accessPoints(~apIdxs);
            obj.receiverAps = accessPoints(apIdxs);
            
            % TODO: solve correctly, not just numel(accessPoints)???
            rcs = sqrt(rcsVar/2)*(randn(numel(accessPoints),1)+1i*randn(numel(accessPoints),1));
            obj.target = Target(targetLoc, rcs);
        end

        function [] = plotSystem(obj, figureHandle)
            if nargin < 2
                figureHandle = figure;
            end

            % TODO:
            % lös med nya properties receiver ap grejen!
            tApsLocs = [obj.transmitAps.pos];
            receiverLocs = [obj.receiverAps.pos];
            usersLocs = [obj.users([1:(obj.adversaryIdx-1), (obj.adversaryIdx+1):end]).pos];
            %obj.transmitAps = accessPoints([1:(receiverIdx-1), (receiverIdx+1):end]);
            %allUsers
            adversaryLoc = obj.users(obj.adversaryIdx).pos;
            targetLoc = obj.target.pos;
            aggrLoc = [tApsLocs, receiverLocs, usersLocs, targetLoc];
            maxX = 1.5*max(aggrLoc(1,:));
            minX = 1.5*min(aggrLoc(1,:));
            maxY = 1.5*max(aggrLoc(2,:));
            minY = 1.5*min(aggrLoc(2,:));

            figure(figureHandle);
            hold on
    
            for i = 1:width(tApsLocs)
                p1 = plot(tApsLocs(1,i), tApsLocs(2,i), 'bx', 'LineWidth', 3, 'MarkerSize', 10);
                %text(tApsLocs(1,i), tApsLocs(2,i), num2str(i), "FontSize", 20);
            end
            
            p2 = plot(receiverLocs(1,:), receiverLocs(2,:), 'bsquare', 'LineWidth', 2, 'MarkerSize', 10);

            p3 = plot(usersLocs(1,:), usersLocs(2,:), 'r*', 'LineWidth', 3, 'MarkerSize', 10);

            p4 = plot(adversaryLoc(1), adversaryLoc(2), 'rdiamond', 'LineWidth', 2, 'MarkerSize', 10);

            p5 = plot(targetLoc(1), targetLoc(2), 'gpentagram', 'LineWidth', 2, 'MarkerSize', 15);
            
            legend([p1, p2, p3, p4, p5], ["transmitter APs", "receiver AP", "User", "Adversary", "Target"]);
            %axis([minX maxX minY maxY])
            axis([-500 500 -500 500])
            hold off
            set(gca,'XTick',[])
            set(gca,'YTick',[])
        end

        function [SINRs] = computeSensingSINR(obj)
            signal = [obj.cSignal;
                      obj.sSignal];



            M = obj.transmitAps.nAntennas;
            nTAps = numel(obj.transmitAps);
            nRAps = numel(obj.receiverAps);

            ATot = zeros(M*(nTAps));

            ATotTot = zeros([size(ATot), nRAps]);
            for k = 1:nRAps
                for i = 1:nTAps
                    for j = 1:nTAps
                        % steering vectors: AP1, AP2, and receiver
                        % plus pi since supposed to be from target to APs according to paper?
                        a1 = obj.transmitAps(i).calcSteeringVector(obj.target, "from");
                        a2 = obj.transmitAps(j).calcSteeringVector(obj.target, "from");
                        ar = obj.receiverAps(k).calcSteeringVector(obj.target, "from");
            
                        A = sqrt(obj.sBetas(i, k)*obj.sBetas(j, k))*...
                            conj(a1)*ar'*(obj.target.rcs(i)*conj(obj.target.rcs(j)))*ar*a2.';
                        ATot(((i-1)*M+1):(i*M),((j-1)*M+1):(j*M)) = A;
                    end
                end
                ATotTot(:,:,k) = ATot;
            end

            sinrSTot = 0;
            for j = 1:nRAps
                sinrSTot = sinrSTot + trace(real(signal'*obj.precoders'*ATotTot(:,:,j)*obj.precoders*signal)...
                    /(width(signal)*nRAps*M*obj.noiseVar));
            end
            SINRs = sinrSTot;
        end

        % TODO:
        % input validation!
        [] = generateSignal(obj, nSymbols, modulationOrder)
        [] = generateChannelProperties(obj, plExp)
        [] = optimizeConfiguration(obj, maxIts, threshold, maxPower, ueSinrConstraint);
        [] = computePrecoders(obj, maxIts, threshold, maxPower, ueSinrConstraint)
        [] = estimateChannelCovariance(obj, maxIts, threshold, channelEstVar)
        [] = estimateTargetAngle(obj)
        [pos] = estimateTargetPositionGrid(obj, cellSize)
        [pos] = estimateTargetPositionMmse(obj, maxIts, threshold, stepSize, options)
    end
end

