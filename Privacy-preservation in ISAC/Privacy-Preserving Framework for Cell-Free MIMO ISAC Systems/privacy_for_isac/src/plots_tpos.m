clear; clc; close all;

%%

clear;
figure;
hold on


%%% DE HÄR ÄR MED ALLT SLUMPAT (tror jag?)

% baseline
f = "data-new-new/b-tpos-2.mat";

load(f);
parameter = distances;
%f = "data_new/uesinr_14.mat";

[pd, sensingSinr] = computeStats(f);
yyaxis left
p1 = plot(parameter, pd, '-x', LineWidth=1.5, MarkerSize=10);
ylim([0,1]);
ylabel("Probability of detection, $P_D$", Interpreter="latex");

yyaxis right
p2 = plot(parameter, 10*log10(sensingSinr), '--x', LineWidth=1.5, MarkerSize=10);
ylim([0, 20]);
ylabel("Sensinr SINR", Interpreter="latex");


% old
f = "data-new-new/o-tpos-2.mat";

load(f);
parameter = distances;
%f = "data_new/uesinr_14.mat";

[pd, sensingSinr] = computeStats(f);
yyaxis left
p1 = plot(parameter, pd, '-*', LineWidth=1.5, MarkerSize=10);
ylim([0,1]);
ylabel("Probability of detection, $P_D$", Interpreter="latex");

yyaxis right
p2 = plot(parameter, 10*log10(sensingSinr), '--*', LineWidth=1.5, MarkerSize=10);
ylim([0, 20]);
ylabel("Sensinr SINR", Interpreter="latex");


% new
f = "data-new-new/n-tpos-2.mat";

load(f);
parameter = distances;
%f = "data_new/uesinr_14_op.mat";

[pd, sensingSinr] = computeStats(f);
yyaxis left
p3 = plot(parameter, pd, '-+', LineWidth=1.5, MarkerSize=10);
ylim([0,1]);
ylabel("Probability of detection, $P_D$", Interpreter="latex");

yyaxis right
p4 = plot(parameter, 10*log10(sensingSinr), '--+', LineWidth=1.5, MarkerSize=10);
ylim([0, 20]);
ylabel("Sensinr SINR", Interpreter="latex");

legend(["baseline", "old", "new", "baseline", "old", "new"]);

xlabel("Distance between adversary and target $[m]$", Interpreter="latex");


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