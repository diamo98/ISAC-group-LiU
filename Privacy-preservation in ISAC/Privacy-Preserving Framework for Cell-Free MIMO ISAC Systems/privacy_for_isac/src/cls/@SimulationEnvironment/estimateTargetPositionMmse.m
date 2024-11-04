function [pos] = estimateTargetPositionMmse(obj, maxIts, threshold, stepSize, options)
    arguments
        obj
        maxIts
        threshold
        stepSize
        options.plot (1, 1) logical = 0
    end
%ESTIMATETARGETPOSITIONMMSE Summary of this function goes here
%   Detailed explanation goes here

% TODO: mycket duplicerad kod från andra targetPosition estimate funktionen
% eftersom linjerna måste skapas av båda algoritmerna så kanske skapa en
% hjälpfunktion för linjer?

    %% parameters
    areaSize = 500;% single side area size so double +- areaSize in each direction
    kMax = 50;
    resolution = 10000;

    linesMat = NaN*ones(resolution, numel(obj.transmitAps)*2);

    e = stepSize; % step size
    b = 10; % batch size, chosen batches are estimate +-e/(b/2) with b/2 steps in each direction,
    % i.e b X b surface with e/(b/2) between each point with estimate point
    % (origo of the gradient batch), at the middle of the b X b matrix

    %% setup lines
    for i = 1:numel(obj.transmitAps)
        x = -areaSize+((areaSize/(resolution/2))):(areaSize/(resolution/2)):areaSize;
        angle = obj.estimatedAngles(i);

        apPos = obj.transmitAps(i).pos;

        k = tan(angle);
        if abs(k) > kMax
            y = x;
            k = tan(angle + pi/2);
            m = apPos(1) + k*apPos(2);
            x = -k*y + m;
        else
            m = apPos(2) - k*apPos(1);
            y = k*x + m;
        end
 
        % valid idx not working properly, case when wrong angle, see
        % varying_h_var_random_system.mat for erronuous case
        % if apPos(2) == obj.target.pos(2)
        %     if apPos(1) < obj.target.pos(1)
        %         validIdx = find(abs(x) <= areaSize & x > apPos(1));
        %     else
        %         validIdx = find(abs(x) <= areaSize & x < apPos(1));
        %     end
        % elseif apPos(2) < obj.target.pos(2)
        %     validIdx = find(abs(y) <= areaSize & y > apPos(2));
        % else
        %     validIdx = find(abs(y) <= areaSize & y < apPos(2));
        % end

        % linesMat(validIdx, ((i*2)-1):(i*2)) = [x(validIdx)', y(validIdx)'];
        % y = y(validIdx);
        % x = x(validIdx);
        linesMat(:, ((i*2)-1):(i*2)) = [x', y'];
    end


    %% Gradient descent
    estimate = [randi([-areaSize, areaSize]), randi([-areaSize, areaSize])]; % init to random


    plotPoints = [estimate];

    
    for i = 1:maxIts
        prevEstimate = estimate;
        [batchX, batchY] = meshgrid(estimate(1)-e+e/(b/2):e/(b/2):estimate(1)+e, estimate(2)-e+e/(b/2):e/(b/2):estimate(2)+e);
        mse = zeros(b);
        
        for j = 1:b^2
            mse(j) = meanSquaredError([batchX(j), batchY(j)], linesMat);
        end

        [FX, FY] = gradient(mse, e/(b/2));
        
        dir = [FX((b^2)/2), FY((b^2)/2)];%./norm([FX((b^2)/2), FY((b^2)/2)]);

        estimate = estimate - e*dir;
        % norm(estimate - prevEstimate)
        if norm(estimate - prevEstimate) < threshold
            break
        end
        plotPoints = [plotPoints; estimate];
    end
    pos = estimate;

    %% plots
    if options.plot
        figure
        hold on;
    
        for i = 1:(width(linesMat)/2)
            plot(linesMat(:,((i*2)-1)), linesMat(:,(i*2)));
        end
        plot(plotPoints(:,1), plotPoints(:,2), 'b+', 'LineWidth', 3);
        plot(plotPoints(:,1), plotPoints(:,2), 'b--', 'LineWidth', 2);
    
        searchR = 40;
        xCirc = searchR * cos((0:resolution).*(2*pi/resolution)) + estimate(1);
        yCirc = searchR * sin((0:resolution).*(2*pi/resolution)) + estimate(2);
        plot(xCirc, yCirc, 'LineWidth', 2);
    
        hold off;
        axis([-areaSize, areaSize, -areaSize, areaSize]);
    end
end

