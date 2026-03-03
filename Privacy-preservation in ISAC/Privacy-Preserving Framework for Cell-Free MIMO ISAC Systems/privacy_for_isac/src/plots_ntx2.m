clear; clc; close all;

%%

clear;
figure(1);
hold on
figure(2);
hold on


%%% DE HÄR ÄR MED ALLT SLUMPAT (tror jag?)
% baseline
%f = "../new-precoder-mc-run1.mat";
f = "data-new-new/b-ntx-1.mat";

load(f);
[pd, sensingSinr] = computeStats(f);
figure(1)
p3 = plot(nrTxApsList, pd, '-.', LineWidth=3, MarkerSize=10);

figure(2)
p4 = plot(nrTxApsList, 10*log10(sensingSinr), '-.', LineWidth=3, MarkerSize=10);

% old
%f = "old-precoder-selection-mc-run1.mat";
f = "data-new-new/o-ntx-1.mat";

load(f);
[pd, sensingSinr] = computeStats(f);
figure(1)
p1 = plot(nrTxApsList, pd, '-*', LineWidth=3, MarkerSize=10);

figure(2)
p2 = plot(nrTxApsList, 10*log10(sensingSinr), '-*', LineWidth=3, MarkerSize=10);

% new
%f = "../new-precoder-mc-run1.mat";
f = "data-new-new/n-ntx-1.mat";

load(f);
[pd, sensingSinr] = computeStats(f);
figure(1)
p3 = plot(nrTxApsList, pd, '-+', LineWidth=3, MarkerSize=10);

figure(2)
p4 = plot(nrTxApsList, 10*log10(sensingSinr), '-+', LineWidth=3, MarkerSize=10);

% combined
%f = "../new-precoder-mc-run1.mat";
f = "data-new-new/c-ntx-1.mat";

load(f);
[pd, sensingSinr] = computeStats(f);
figure(1)
p3 = plot(nrTxApsList, pd, '-x', LineWidth=3, MarkerSize=10);

figure(2)
p4 = plot(nrTxApsList, 10*log10(sensingSinr), '-x', LineWidth=3, MarkerSize=10);

figure(1)
ylabel("Probability of detection, $P_D$", Interpreter="latex");
legend(["Naive precoder", "Naive precoder + sel", "MI-precoder", "MI-precoder + sel"])
xlabel("Number of transmitting access points", Interpreter="latex");
xlim([4, 6]);
ax = gca;
ax.XTick = unique(round(ax.XTick));
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15);
hold off;
box on;

figure(2)
ylabel("Sensinr SINR [dB]", Interpreter="latex");
legend(["Naive precoder", "Naive precoder + sel", "MI-precoder", "MI-precoder + sel"])
xlabel("Number of transmitting access points", Interpreter="latex");
xlim([4, 6]);
ax = gca;
ax.XTick = unique(round(ax.XTick));
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15);
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