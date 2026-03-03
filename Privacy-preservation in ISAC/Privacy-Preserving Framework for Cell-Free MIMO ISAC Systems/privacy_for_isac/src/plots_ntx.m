clear; clc; close all;

%%

clear;
figure;
hold on


%%% DE HÄR ÄR MED ALLT SLUMPAT (tror jag?)
% old
%f = "old-precoder-selection-mc-run1.mat";
f = "data-new-new/o-ntx-1.mat";

load(f);
[pd, sensingSinr] = computeStats(f);
yyaxis left
p1 = plot(nrTxApsList, pd, '-*', LineWidth=1.5, MarkerSize=10);
ylim([0,1]);
ylabel("Probability of detection, $P_D$", Interpreter="latex");

yyaxis right
p2 = plot(nrTxApsList, 10*log10(sensingSinr), '--*', LineWidth=1.5, MarkerSize=10);
ylim([0, 20]);
ylabel("Sensinr SINR", Interpreter="latex");


% new
%f = "../new-precoder-mc-run1.mat";
f = "data-new-new/n-ntx-1.mat";

load(f);
[pd, sensingSinr] = computeStats(f);
yyaxis left
p3 = plot(nrTxApsList, pd, '-+', LineWidth=1.5, MarkerSize=10);
ylim([0,1]);
ylabel("Probability of detection, $P_D$", Interpreter="latex");

yyaxis right
p4 = plot(nrTxApsList, 10*log10(sensingSinr), '--+', LineWidth=1.5, MarkerSize=10);
ylim([0, 20]);
ylabel("Sensinr SINR", Interpreter="latex");


% baseline
%f = "../new-precoder-mc-run1.mat";
f = "data-new-new/b-ntx-1.mat";

load(f);
[pd, sensingSinr] = computeStats(f);
yyaxis left
p3 = plot(nrTxApsList, pd, '-.', LineWidth=1.5, MarkerSize=10);
ylim([0,1]);
ylabel("Probability of detection, $P_D$", Interpreter="latex");

yyaxis right
p4 = plot(nrTxApsList, 10*log10(sensingSinr), '--.', LineWidth=1.5, MarkerSize=10);
ylim([0, 20]);
ylabel("Sensinr SINR", Interpreter="latex");


% combined
%f = "../new-precoder-mc-run1.mat";
f = "data-new-new/c-ntx-1.mat";

load(f);
[pd, sensingSinr] = computeStats(f);
yyaxis left
p3 = plot(nrTxApsList, pd, '-x', LineWidth=1.5, MarkerSize=10);
ylim([0,1]);
ylabel("Probability of detection, $P_D$", Interpreter="latex");

yyaxis right
p4 = plot(nrTxApsList, 10*log10(sensingSinr), '--x', LineWidth=1.5, MarkerSize=10);
ylim([0, 30]);
ylabel("Sensinr SINR", Interpreter="latex");

legend(["old PD", "new PD", "baseline PD", "combined PD", "old sSINR", "new sSINR", "baseline sSINR", "combined sSINR"]);

xlabel("Number of transmitting APs", Interpreter="latex");


hold off;
box on;

%%

function [pd, sensingSinr] = computeStats(file)
    load(file)

    r = cellSize/sqrt(pi);
    
    dataPoints = length(mmseResultsCell);
    correctMmse = zeros(1, dataPoints);
    sensingSinr = zeros(1, dataPoints);

    for j = 1:dataPoints
        mmseResults = mmseResultsCell{j};
        targetLocs = targetLocsCell{j};
        sensingSinr(j) = sum(sensingSinrResultsCell{j})/mcRuns;
    
        for i = 1:mcRuns
        
            mmseEstimate = mmseResults(i,:);
            targetLoc = targetLocs(i,:);
        
            if (targetLoc(1) - mmseEstimate(1))^2 + (targetLoc(2) - mmseEstimate(2))^2 <= r^2
                correctMmse(j) = correctMmse(j) + 1;
            end
        end
    end
    pd = correctMmse./mcRuns;
end