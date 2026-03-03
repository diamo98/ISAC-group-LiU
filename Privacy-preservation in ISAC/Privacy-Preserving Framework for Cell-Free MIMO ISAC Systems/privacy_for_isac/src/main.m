clc; clear; close all;

% TODO: input-validering på alla klasser och refaktorera lite typ precoder
% default i AccessPoint eller i SimulationEnvironment när APn:a de skapas?

% TODO: Gör lite structs i SimulationEnvironment med t.ex kanalmatriser
% tillsammans med användarna de avser i enstruct typ
% struct UserData { User u; Mt X Mr complex double channelMatrix; }
% Kanske till och med en matris per mottagar antenn för ifall det skulle
% behövas senare + blir lite snyggare i estimateChannelCovariance

% TODO: Kanske ändra struktur för lite konstigt att ändra på ett tillstånd
% anonymt hela tiden?

%% monte carlo trials random system model
% addpath("../../cvx");
% run("cvx_setup.m");
% addpath("cls");
% addpath("functions");

clear;

%nUsersList = [1,2,3,4,5,6];
%nUsers = 3;
%nAps = 4;
%nReceiversList = 1;%[1, 2, 3];
apsLocs = [-300, -300;
           -300, 300;
           300, -300;
           300, 300];

% apsLocsList = [-300, -300;
%                -300, 300;
%                300, -300;
%                300, 300
%                0, 300;
%                0, -300;
%                300, 0;
%                -300, 0];

% nrTxApsList = [4, 5, 6, 7, 8];

targetLoc = [0, 0];

usersLocs = [100, 150;
             -150, 50;
             0, -100];
% distances = [10, 20, 50, 100];
% targetLocList = [90, 150;
%                  80, 150;
%                  50, 150
%                  0, 150];
%nReceivers = 1;
receiverIdxs = 1;
adversaryIdx = 1;
gridSize = 500;

nSymbols = 16;
%nSymbolsList = 2.^(0:9);
modulationOrder = 16;

maxEmIts = 100;
maxCvxIts = 2; % var 100 men borde räcka med 1? ty cvx. 2 för att få en hyffsad Taylor också
maxMmseIts = 100;

emThreshold = 1e-3;
cvxThreshold = 1;
mmseThreshold = 1;

apAntennas = 64;
%apAntennasList = 2.^(4:8);
userAntennas = 4;

nObservations = 32;


%maxPower = apAntennas * (10^(35/10)/1000);
maxPower = (10^(58/10)/1000);
%maxPower = apAntennas * (10^(40/10)/1000);
%maxPowerList = apAntennas * (10.^((30:5:50)./10)./1000);

%sensingSinr = 10^(11/10);
%ueSinr = 10^(-2/10);
ueSinr = 10^(5/10);
fractionSensingUeSinr = [0.5, 1, 2, 4, 10, 15, 20];

%ueSinrList = 10.^((-2:2)./10);

channelEstVar = (10^(-40/10)/1000); % var -75 dB
%channelEstVarList = 10.^((-80:1:-60)./10);

rcsVar = 10^(10/10);
noiseVar = 10^(-94/100)/1000;

cellSize = 10;
stepSize = 0.75;

mcRuns = 1;

mmseResultsCell = {};
gridResultsCell = {};
anglesResultsCell = {};
sensingSinrResultsCell = {};
targetLocsCell = {};
cvxCvgCell = {};

for fracSinr = fractionSensingUeSinr
    sensingSinr = fracSinr*ueSinr;

    mmseResults = [];
    gridResults = [];
    anglesResults = [];
    sensingSinrResults = [];
    targetLocs = [];
    cvxCvg = [];
    for i = 1:mcRuns
        disp("bokstav")
        sensingSinr
        i
        tic
        %[usersLocs, apsLocs, targetLoc] = randomizeSystem(nUsers, nAps, gridSize);
        %adversaryIdx = randi([1, nUsers]);
        %receiverIdxs = randperm(nAps, nReceivers);
        cvg = false;
        cvg_count = 0;
        cvg_count_max = 100;
        while ~cvg && cvg_count < 100
            simEnv = SimulationEnvironment(apAntennas, userAntennas, ...
                            apsLocs, receiverIdxs, usersLocs, adversaryIdx, targetLoc, ...
                            rcsVar, noiseVar, maxPower, nObservations);
            simEnv.generateSignal(nSymbols, modulationOrder);
            simEnv.generateChannelProperties();
            % simEnv.optimizeConfiguration(maxCvxIts, cvxThreshold, maxPower, ueSinr);
            % cvg = simEnv.computePrecoders(maxCvxIts, cvxThreshold, maxPower, ueSinr);
            cvg = simEnv.computePrecoders2(maxCvxIts, cvxThreshold, maxPower, ueSinr, sensingSinr);
            cvg_count = cvg_count + 1;
        end

        simEnv.estimateChannelCovariance(maxEmIts, emThreshold, channelEstVar);
        simEnv.estimateTargetAngle();
        posEstGrid = simEnv.estimateTargetPositionGrid(cellSize);

        %plot?
        posEstMmse = simEnv.estimateTargetPositionMmse(maxMmseIts, mmseThreshold, stepSize, "plot",false);
        %simEnv.plotSystem(gcf);

        anglesResults = [anglesResults; simEnv.estimatedAngles];
        gridResults = [gridResults; posEstGrid];
        mmseResults = [mmseResults; posEstMmse];
        sensingSinrResults = [sensingSinrResults; simEnv.computeSensingSINR()];
        targetLocs = [targetLocs; targetLoc];
        cvxCvg = [cvxCvg; cvg];
        toc
    end
    anglesResultsCell = [anglesResultsCell, {anglesResults}];
    mmseResultsCell = [mmseResultsCell, {mmseResults}];
    gridResultsCell = [gridResultsCell, {gridResults}];
    sensingSinrResultsCell = [sensingSinrResultsCell, {sensingSinrResults}];
    targetLocsCell = [targetLocsCell, {targetLocs}];
    cvxCvgCell = [cvxCvgCell, {cvxCvg}]; 
end

%save("o-quot-1.mat");
