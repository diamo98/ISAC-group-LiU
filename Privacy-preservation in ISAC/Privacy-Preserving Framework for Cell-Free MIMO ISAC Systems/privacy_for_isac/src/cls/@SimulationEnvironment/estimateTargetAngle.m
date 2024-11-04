function [] = estimateTargetAngle(obj)
%ESTIMATETARGETPOSITION Summary of this function goes here
%   Detailed explanation goes here

% TODO:
% 1: 
% hela funktionen är ganska felplacerad, konstigt att centrala styrenheten estimerar
% kovariansen utifrån en användares perspektiv, plus världens mest
% temporära ickevälskrivna fuskfunktion
% 2:
% ta in grid och kanske AP-vinklar som argument

angleSteps = 1000;

%B = zeros([2*numel(obj.transmitAps), angleSteps]);
%cB = zeros([2*numel(obj.transmitAps), angleSteps]);

totSig = [obj.cSignal;
          obj.sSignal];
nSymbols = width(totSig);

for i = 1:numel(obj.transmitAps)
    ap = obj.transmitAps(i);
    correctX = ap.precoderMatrix*totSig;

    % angle sweep dependent on what rotated quadrant the ap is in, ap at 0,0 wont 
    % provide any useful results anyways.
    % Aps at the lines y=x and y=-x wont matter which of its two available ifs they fall
    % under since the result should be the same anyways.
    % Since cos in the exponent of array response vector, its only unique
    % values are 0-180 so avoid crossing over those to get a unique result

    % only sensible way to choose angle since periodicity of array response
    % vector cosine 0-180 only unique values
    if ap.pos(2) >= obj.target.pos(2)
        angleSweep = linspace(0, -pi, angleSteps);
    else
        angleSweep = linspace(0, pi, angleSteps);
    end

    tmpCor = [];
    tmpEst = [];

    R = correctX*correctX';

    meanEstR = 0;
    for k = 1:nSymbols
        % normalized
        meanEstR = meanEstR + obj.estimatedR{i}{k}./max(obj.estimatedR{i}{k}(:));
    end
    meanEstR = (1/nSymbols)*meanEstR;

    for j = 1:angleSteps
        a = ap.calcSteeringVector(angleSweep(j));

        tmpEst = [tmpEst; real(a'*meanEstR*a)];
        tmpCor = [tmpCor; real(a'*R*a)];
    end
    
%     meanAngle = 0;
%     for k = 1:nSymbols
%         [~, idx] = max(tmpEst(:, k));
%         
% %         figure
% %         plot(tmpEst(:,k));
% 
%         meanAngle = meanAngle + angleSweep(idx);
%     end
%     obj.estimatedAngles(i) = meanAngle/nSymbols;
        %B(2*i,:) = angleSweep;
        %cB(2*i,:) = angleSweep;
    [~, idx] = max(tmpEst);
    obj.estimatedAngles(i) = angleSweep(idx);
    [~, idx] = max(tmpCor);
    obj.correctAngles(i) = angleSweep(idx);

    %{
    figure
    plot(angleSweep.*(180/pi), tmpCor, 'LineWidth', 2);
    xlabel("Angle in Degrees")
    set(findall(gcf,'-property','FontSize'),'FontSize',13);
    title("correct");

    figure
    plot(angleSweep.*(180/pi), tmpEst, 'LineWidth', 2);
    xlabel("Angle in Degrees")
    set(findall(gcf,'-property','FontSize'),'FontSize',13);
    title("estimate");
    %}
end
%obj.correctAngles
%obj.estimatedAngles

% TODO:
% använd aBa funktionen och svep över t.ex -90, 90 grader (utifrån vart de
% är positionerade gentemot målet?) gör samma sak med den egentliga
% precodern och se om de matchar eller om vinkeln är åt rätt håll.

% gör rutnät kanske -500, 500 och följ vinklarna och se om de träffar någon
% gemensam ruta.

% de plots som är av intresse är sensing SINR mot SINR_ui constrainten


end

