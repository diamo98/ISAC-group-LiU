clear; clc; close all;

%%
[pd, sensingSinr] = computeStats("/Users/henrikakesson/OneDrive - Linköpings universitet/research/c-quot-2.mat");
sinr = 10*log10(sensingSinr);
%%
clear;
figure(1);
hold on
% figure(2);
% hold on

folder = "C:\Users\henak56\research-data\privacy_cell_free_isac\old\data4\";
% baseline (utan något förbättrande, sekretessmässigt)
f = append(folder, "b-quot-2.mat");

load(f);
[pd, sensingSinr] = computeStats(f);
if isscalar(pd)
    pd = pd*ones(1, numel(fractionSensingUeSinr));
end
if isscalar(sensingSinr)
    sensingSinr = sensingSinr*ones(1, numel(fractionSensingUeSinr));
end
figure(1)
yyaxis left
p3 = plot(fractionSensingUeSinr, pd, '-.', LineWidth=3, MarkerSize=10);

figure(1)
yyaxis right
p4 = plot(fractionSensingUeSinr, 10*log10(sensingSinr), '-.', LineWidth=3, MarkerSize=10);

% old (gamla precodern, med selection)
f = append(folder, "o-quot-2.mat");

load(f);
[pd, sensingSinr] = computeStats(f);
if isscalar(pd)
    pd = pd*ones(1, numel(fractionSensingUeSinr));
end
if isscalar(sensingSinr)
    sensingSinr = sensingSinr*ones(1, numel(fractionSensingUeSinr));
end
figure(1)
yyaxis left
p1 = plot(fractionSensingUeSinr, pd, '-*', LineWidth=3, MarkerSize=10);

figure(1)
yyaxis right
p2 = plot(fractionSensingUeSinr, 10*log10(sensingSinr), '-*', LineWidth=3, MarkerSize=10);

% new (ny precoder, utan selection)
f = append(folder, "n-quot-2.mat");

load(f);
[pd, sensingSinr] = computeStats(f);
figure(1)
yyaxis left
p3 = plot(fractionSensingUeSinr, pd, '-+', LineWidth=3, MarkerSize=10);

figure(1)
yyaxis right
p4 = plot(fractionSensingUeSinr, 10*log10(sensingSinr), '-+', LineWidth=3, MarkerSize=10);

% combined (ny precoder, med selection)
f = append(folder, "c-quot-1.mat");

load(f);
[pd, sensingSinr] = computeStats(f);
figure(1)
yyaxis left
p3 = plot(fractionSensingUeSinr, pd, '-x', LineWidth=3, MarkerSize=10);

figure(1)
yyaxis right
p4 = plot(fractionSensingUeSinr, 10*log10(sensingSinr), '-x', LineWidth=3, MarkerSize=10);

figure(1)
yyaxis left
ylabel("Probability of detection, $P_D$", Interpreter="latex");
legend(["b pd", "o pd", "n pd", "c pd", "b sinr", "o sinr", "n sinr", "c sinr"])
xlabel("SINR constraint fraction $\frac{\gamma_s}{\gamma_{UE}}$", Interpreter="latex");
xlim([0, 15]);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15);
yyaxis right
ylabel("Sensinr SINR [dB]", Interpreter="latex")
hold off;
box on;

% figure(2)
% ylabel("Sensinr SINR [dB]", Interpreter="latex");
% legend(["Naive precoder", "Naive precoder + sel", "MI-precoder", "MI-precoder + sel"])
% xlabel("SINR constraint fraction $\frac{\gamma_s}{\gamma_{UE}}$", Interpreter="latex");
% xlim([0, 20]);
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15);
% hold off;
% box on;

matlab2tikz('filename', 'tikfigure.tex', 'figurehandle', figure(1));

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
        cvg = (cvxCvgCell{j} == 1);
        converged_runs = sum(cvg);
        sensingSinr(j) = sum(sensingSinrResultsCell{j}(cvg))/converged_runs;
        for i = 1:mcRuns
            if cvg(i)
                mmseEstimate = mmseResults(i,:);
                targetLoc = targetLocs(i,:);
            
                if (targetLoc(1) - mmseEstimate(1))^2 + (targetLoc(2) - mmseEstimate(2))^2 <= r^2
                    correctMmse(j) = correctMmse(j) + 1;
                end
            end
        end
    end
    pd = correctMmse./converged_runs;
end