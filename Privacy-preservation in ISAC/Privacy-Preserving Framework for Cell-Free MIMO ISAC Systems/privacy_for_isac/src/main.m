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
clear;

%nUsersList = [1,2,3,4,5,6];
nUsers = 3;
nAps = 6;
nReceiversList = 4;%[1, 2, 3];
%nReceivers = 3;
gridSize = 500;

nSymbols = 16;
%nSymbolsList = 2.^(0:9);
modulationOrder = 16;

maxEmIts = 100;
maxCvxIts = 100;
maxMmseIts = 100;

emThreshold = 1e-3;
cvxThreshold = 1;
mmseThreshold = 1;

apAntennas = 64;
%apAntennasList = 2.^(4:8);
userAntennas = 4;

nObservations = 32;

maxPower = apAntennas * (10^(35/10)/1000);
%maxPower = apAntennas * (10^(40/10)/1000);
%maxPowerList = apAntennas * (10.^((30:5:50)./10)./1000);

ueSinr = 10^(2/10);
%ueSinrList = 10.^((0:5)./10);


channelEstVar = 10^(-75/10);
%channelEstVarList = 10.^((-80:1:-60)./10);

rcsVar = 10^(10/10);
noiseVar = 10^(-94/100)/1000;

cellSize = 10;
stepSize = 1;

mcRuns = 1;

mmseResultsCell = {};
gridResultsCell = {};
anglesResultsCell = {};
sensingSinrResultsCell = {};
targetLocsCell = {};

for nReceivers = nReceiversList
    mmseResults = [];
    gridResults = [];
    anglesResults = [];
    sensingSinrResults = [];
    targetLocs = [];
    for i = 1:mcRuns
        nReceivers
        i
        tic
        [usersLocs, apsLocs, targetLoc] = randomizeSystem(nUsers, nAps, gridSize);
        adversaryIdx = randi([1, nUsers]);
        receiverIdxs = randperm(nAps, nReceivers);

        simEnv = SimulationEnvironment(apAntennas, userAntennas, ...
                        apsLocs, receiverIdxs, usersLocs, adversaryIdx, targetLoc, ...
                        rcsVar, noiseVar, maxPower, nObservations);
        simEnv.generateSignal(nSymbols, modulationOrder);
        simEnv.generateChannelProperties();
        simEnv.optimizeConfiguration(maxCvxIts, cvxThreshold, maxPower, ueSinr);
        simEnv.computePrecoders(maxCvxIts, cvxThreshold, maxPower, ueSinr);
        simEnv.estimateChannelCovariance(maxEmIts, emThreshold, channelEstVar);
        simEnv.estimateTargetAngle();
        posEstGrid = simEnv.estimateTargetPositionGrid(cellSize);
        posEstMmse = simEnv.estimateTargetPositionMmse(maxMmseIts, mmseThreshold, stepSize);
        
        anglesResults = [anglesResults; simEnv.estimatedAngles];
        gridResults = [gridResults; posEstGrid];
        mmseResults = [mmseResults; posEstMmse];
        sensingSinrResults = [sensingSinrResults; simEnv.computeSensingSINR()];
        targetLocs = [targetLocs; targetLoc];
        toc
    end
    anglesResultsCell = [anglesResultsCell, {anglesResults}];
    mmseResultsCell = [mmseResultsCell, {mmseResults}];
    gridResultsCell = [gridResultsCell, {gridResults}];
    sensingSinrResultsCell = [sensingSinrResultsCell, {sensingSinrResults}];
    targetLocsCell = [targetLocsCell, {targetLocs}];
end

save("data_new/nraps_16_op_single.mat");

%% P_D plots
clear;
figure;
hold on

load("data_new/nraps_14.mat");

r = cellSize/sqrt(pi);

dataPoints = length(mmseResultsCell);
correctMmse = zeros(1, dataPoints);
correctGrid = zeros(1, dataPoints);
sensingSinr1 = zeros(1, dataPoints);

for j = 1:dataPoints
    mmseResults = mmseResultsCell{j};
    gridResults = gridResultsCell{j};
    targetLocs = targetLocsCell{j};
    % sensingSinr1(j) = sum(sensingSinrResultsCell{j})/mcRuns;

    for i = 1:mcRuns
    
        mmseEstimate = mmseResults(i,:);
        gridEstimate = gridResults(i,:);
        targetLoc = targetLocs(i,:);
    
        if (targetLoc(1) - mmseEstimate(1))^2 + (targetLoc(2) - mmseEstimate(2))^2 <= r^2
            correctMmse(j) = correctMmse(j) + 1;
        end
    
        % inte helt rättvist för då är alltid gissningen i mitten på celler
        % och då kan två gissnigar ge rätt
        if (all(abs(targetLoc - gridEstimate) <= cellSize/2))
            correctGrid(j) = correctGrid(j) + 1;
        end
    end
end

% plot(10*log10(maxPowerList*1000/apAntennas), correctMmse./mcRuns, '--+', LineWidth=1.5);
plot(nReceiversList, correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(nUsersList, correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(10*log10(ueSinrList), correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(10*log10(channelEstVarList), correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(nSymbolsList, correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(apAntennasList, correctMmse./mcRuns, '--+', LineWidth=1.5);

%TRADEOFF PLOTS
% plot(10*log10(sensingSinr1), correctMmse./mcRuns, '--*', LineWidth=1.5, MarkerSize=10);
% point_labels = [];
% y_offset = 0.05;
% x_offset = -0.5;
% for i = 1:dataPoints
%     point_labels = [point_labels, sprintf("N_{Rx}=%d", i)];
% end
% text(10*log10(sensingSinr1) + x_offset, correctMmse./mcRuns + y_offset, point_labels);


%%% SECOND PLOT

clearvars -except sensingSinr1;
load("data_new/nraps_14_op.mat");

r = cellSize/sqrt(pi);

dataPoints = length(mmseResultsCell);
correctMmse = zeros(1, dataPoints);
correctGrid = zeros(1, dataPoints);
sensingSinr2 = zeros(1, dataPoints);

for j = 1:dataPoints
    mmseResults = mmseResultsCell{j};
    gridResults = gridResultsCell{j};
    targetLocs = targetLocsCell{j};
    % sensingSinr2(j) = sum(sensingSinrResultsCell{j})/mcRuns;

    for i = 1:mcRuns
    
        mmseEstimate = mmseResults(i,:);
        gridEstimate = gridResults(i,:);
        targetLoc = targetLocs(i,:);
    
        if (targetLoc(1) - mmseEstimate(1))^2 + (targetLoc(2) - mmseEstimate(2))^2 <= r^2
            correctMmse(j) = correctMmse(j) + 1;
        end
    
        % inte helt rättvist för då är alltid gissningen i mitten på celler
        % och då kan två gissnigar ge rätt
        if (all(abs(targetLoc - gridEstimate) <= cellSize/2))
            correctGrid(j) = correctGrid(j) + 1;
        end
    end
end

% plot(10*log10(maxPowerList*1000/apAntennas), correctMmse./mcRuns, '--*', LineWidth=1.5);
plot(nReceiversList, correctMmse./mcRuns, '--*', LineWidth=1.5);
% plot(nUsersList, correctMmse./mcRuns, '--*', LineWidth=1.5);
% plot(10*log10(ueSinrList), correctMmse./mcRuns, '--*', LineWidth=1.5);
% plot(10*log10(channelEstVarList), correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(nSymbolsList, correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(apAntennasList, correctMmse./mcRuns, '--+', LineWidth=1.5);

%TRADEOFF PLOTS
% plot(10*log10(sensingSinr2), correctMmse./mcRuns, '--*', LineWidth=1.5, MarkerSize=10);

hold off;
% axis([min(10*log10(maxPowerList*1000/apAntennas)), max(10*log10(maxPowerList*1000/apAntennas)), 0, 1e0]);
% xlabel("Maximum power per antenna element [dBm]");
% ylabel("Probability of detection");

axis([min(nReceiversList), max(nReceiversList), 0, 1e0]);
xlabel("Number of receiver APs");
ylabel("Probability of detection"); 

% axis([min(nUsersList), max(nUsersList), 0, 1e0]);
% xlabel("Number of users");
% ylabel("Probability of detection");

% axis([min(10*log10(ueSinrList)), max(10*log10(ueSinrList)), 0, 1e0]);
% xlabel("Minimum user SINR [dB]");
% ylabel("Probability of detection");

% axis([min(10*log10(channelEstVarList)), max(10*log10(channelEstVarList)), 0, 1e0]);
% xlabel("Channel estimation variance");
% ylabel("Probability of detection");

% axis([min(nSymbolsList), max(nSymbolsList), 0, 1e0]);
% xlabel("Number of symbols");
% ylabel("Probability of detection");

% axis([min(apAntennasList), max(apAntennasList), 0, 1e0]);
% xlabel("Number of AP antennas");
% ylabel("Probability of detection");

legend(["Without AP selection", "With AP selection"]);

%TRADEOFF PLOTS

% axis([min(10*log10([sensingSinr1, sensingSinr2])), max(10*log10([sensingSinr1, sensingSinr2])), 0, 1e0]);
% xlabel("Sensing SINR [dB]")
% ylabel("Probability of detection");
