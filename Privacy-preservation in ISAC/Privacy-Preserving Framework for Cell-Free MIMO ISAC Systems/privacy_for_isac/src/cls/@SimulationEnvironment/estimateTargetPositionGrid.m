function [pos] = estimateTargetPositionGrid(obj, cellSize)
%ESTIMATETARGETPOSITION Summary of this function goes here
%   Detailed explanation goes here
    areaSideSize = 1000; % måste vara delbart med 2 för tillfället
    stepSize = cellSize; % var 50 förut men det är lite väl stort?
    lineRes = 0.1;
    % hold on;
    crossedCells = [];
    angleSweep = linspace(0, pi, 1000);

    for i = 1:numel(obj.transmitAps)
        x = -areaSideSize/2:lineRes:areaSideSize/2;

        angle = obj.estimatedAngles(i);

        apPos = obj.transmitAps(i).pos;

        k = tan(angle);
        if abs(k) > stepSize
            y = x;
            k = tan(angle + pi/2);
            m = apPos(1) + k*apPos(2);
            x = -k*y + m;
            % rotatePoints = true;
        else
            m = apPos(2) - k*apPos(1);
            y = k*x + m;
        end
 
        % valid idx not working properly, case when wrong angle, see
        % varying_h_var_random_system.mat for erronuous case
        % if apPos(2) == obj.target.pos(2)
        %     if apPos(1) < obj.target.pos(1)
        %         validIdx = find(abs(x) <= areaSideSize/2 & x > apPos(1));
        %     else
        %         validIdx = find(abs(x) <= areaSideSize/2 & x < apPos(1));
        %     end
        % elseif apPos(2) < obj.target.pos(2)
        %     validIdx = find(abs(y) <= areaSideSize/2 & y > apPos(2));
        % else
        %     validIdx = find(abs(y) <= areaSideSize/2 & y < apPos(2));
        % end
        % 
        % y = y(validIdx);
        % x = x(validIdx);

        

        % plot(x, y);
        points = unique(round([x',y']./stepSize), 'rows').*stepSize;

        crossedCells = [crossedCells; points];
    end
    % hold off;
    % legend(["1", "2", "3", "4"]);
    % axis([-areaSideSize/2, areaSideSize/2, -areaSideSize/2, areaSideSize/2]);
    count = -inf;
    mostCrossedCell = [nan, nan];

    uniqueCrossedCells = unique(crossedCells, 'rows');
    for i = 1:height(uniqueCrossedCells)

        % vector occurance check
        tmpCount = sum(sum(crossedCells == uniqueCrossedCells(i,:), 2) == 2);
        % TODO: ändra så den slumpar om två har samma?
        if tmpCount > count
            curPoint = uniqueCrossedCells(i,:);
            if curPoint(1) <= areaSideSize/2 && curPoint(1) >= -areaSideSize/2 ...
                && curPoint(2) <= areaSideSize/2 && curPoint(2) >= -areaSideSize/2
                count = tmpCount;
                mostCrossedCell = curPoint;
            end
        end
    end
    pos = mostCrossedCell;
end
