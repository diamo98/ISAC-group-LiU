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
addpath("../../cvx");
run("cvx_setup.m");
addpath("cls");
addpath("functions");

clear;

%nUsersList = [1,2,3,4,5,6];
%nUsers = 3;
%nAps = 4;
%nReceiversList = 1;%[1, 2, 3];
apsLocs = [-300, -300;
           -300, 300;
           300, -300;
           300, 300];
usersLocs = [100, 150;
             -150, 50;
             0, -100];
targetLoc = [0, 0];
%nReceivers = 1;
receiverIdxs = 1;
adversaryIdx = 1;
gridSize = 500;

nSymbols = 16;
%nSymbolsList = 2.^(0:9);
modulationOrder = 16;

maxEmIts = 100;
maxCvxIts = 1; % var 100 men borde räcka med 1? ty cvx
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

%ueSinr = 10^(2/10);
ueSinrList = 10.^((-2:2)./10);

channelEstVar = (10^(-40/10)/1000); % var -75 dB
%channelEstVarList = 10.^((-80:1:-60)./10);

rcsVar = 10^(10/10);
noiseVar = 10^(-94/100)/1000;

cellSize = 10;
stepSize = 0.75;


mcRuns = 100;

mmseResultsCell = {};
gridResultsCell = {};
anglesResultsCell = {};
sensingSinrResultsCell = {};
targetLocsCell = {};
cvxCvgCell = {};

for ueSinr = ueSinrList
    mmseResults = [];
    gridResults = [];
    anglesResults = [];
    sensingSinrResults = [];
    targetLocs = [];
    cvxCvg = [];
    for i = 1:mcRuns
        ueSinr
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
            simEnv.optimizeConfiguration(maxCvxIts, cvxThreshold, maxPower, ueSinr);
            cvg = simEnv.computePrecoders(maxCvxIts, cvxThreshold, maxPower, ueSinr);
            % cvg = simEnv.computePrecoders2(maxCvxIts, cvxThreshold, maxPower, ueSinr, 10^(13/10));
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


save("old-precoder-selection-mc-run1.mat");

%% P_D plots
clear;
figure;
hold on

load("data_new/nraps_18.mat");

r = cellSize/sqrt(pi);

dataPoints = length(mmseResultsCell);
correctMmse = zeros(1, dataPoints);
correctGrid = zeros(1, dataPoints);
sensingSinr = zeros(1, dataPoints);

for j = 1:dataPoints
    mmseResults = mmseResultsCell{j};
    gridResults = gridResultsCell{j};
    targetLocs = targetLocsCell{j};
    sensingSinr(j) = sum(sensingSinrResultsCell{j})/mcRuns;

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

% plot(10*log10(maxPowerList*1000/apAntennas), correctMmse./mcRuns, '-+', LineWidth=1.5, MarkerSize=10);
% plot(nReceiversList, correctMmse./mcRuns, '-+', LineWidth=1.5, MarkerSize=10);
% plot(nApsList, correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(nUsersList, correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(10*log10(ueSinrList), correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(10*log10(channelEstVarList), correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(nSymbolsList, correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(apAntennasList, correctMmse./mcRuns, '--+', LineWidth=1.5);

%TRADEOFF PLOTS
yyaxis left
p1 = plot(nReceiversList, correctMmse./mcRuns, '-*', LineWidth=1.5, MarkerSize=10);
yyaxis right
p2 = plot(nReceiversList, 10*log10(sensingSinr), '--*', LineWidth=1.5, MarkerSize=10);
%%% SECOND PLOT

clearvars -except p1 p2;
load("data_new/nraps_18_op.mat");

r = cellSize/sqrt(pi);

dataPoints = length(mmseResultsCell);
correctMmse = zeros(1, dataPoints);
correctGrid = zeros(1, dataPoints);
sensingSinr = zeros(1, dataPoints);

for j = 1:dataPoints
    mmseResults = mmseResultsCell{j};
    gridResults = gridResultsCell{j};
    targetLocs = targetLocsCell{j};
    sensingSinr(j) = sum(sensingSinrResultsCell{j})/mcRuns;

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

% plot(10*log10(maxPowerList*1000/apAntennas), correctMmse./mcRuns, '--*', LineWidth=1.5, MarkerSize=10);
% plot(nReceiversList, correctMmse./mcRuns, '--*', LineWidth=1.5, MarkerSize=10);
% plot(nApsList, correctMmse./mcRuns, '--*', LineWidth=1.5);
% plot(nUsersList, correctMmse./mcRuns, '--*', LineWidth=1.5);
% plot(10*log10(ueSinrList), correctMmse./mcRuns, '--*', LineWidth=1.5);
% plot(10*log10(channelEstVarList), correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(nSymbolsList, correctMmse./mcRuns, '--+', LineWidth=1.5);
% plot(apAntennasList, correctMmse./mcRuns, '--+', LineWidth=1.5);

%TRADEOFF PLOTS
yyaxis left
p3 = plot(nReceiversList, correctMmse./mcRuns, ':+', LineWidth=1.5, MarkerSize=10);
ylabel("Probability of detection, P_D");
axis([1,4,0,1]);
yyaxis right
p4 = plot(nReceiversList, 10*log10(sensingSinr), '-.+', LineWidth=1.5, MarkerSize=10);
ylabel("Sensing SINR [dB]");
axis([1,3,0,25]);
% legend([p1, p2 ,p3, p4], ...
%     ["No selection: P_d", "No Selection: SINR", "With selection: P_d", "With Selection: SINR"])
legend([p2, p4 ,p1, p3], ...
    ["No Selection: SINR", "With Selection: SINR", "No selection: P_D","With selection: P_D"])

curtick = get(gca, 'XTick');
set(gca, 'XTick', unique(round(curtick)));

hold off;

box on;

% axis([min(10*log10(maxPowerList*1000/apAntennas)), max(10*log10(maxPowerList*1000/apAntennas)), 0, 1e0]);
% xlabel("Maximum power per antenna element [dBm]");
% ylabel("Probability of detection, P_D");

% axis([min(nReceiversList), max(nReceiversList), 0, 1e0]);
% axis([1, 3, 0, 1e0]);
% xlabel("Number of receiver APs");
% ylabel("Probability of detection, P_D");
% curtick = get(gca, 'XTick');
% set(gca, 'XTick', unique(round(curtick)));

% axis([min(nApsList), max(nApsList), 0, 1e0]);
% xlabel("Number of receiver APs");
% ylabel("Probability of detection"); 

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

% legend(["Without AP selection", "With AP selection"]);

%TRADEOFF PLOTS
% xlabel("Number of receiver APs")
