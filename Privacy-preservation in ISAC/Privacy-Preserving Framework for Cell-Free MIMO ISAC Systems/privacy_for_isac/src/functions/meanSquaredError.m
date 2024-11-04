function [mse] = meanSquaredError(point, lines)
%MEANSQUAREDERROR Summary of this function goes here
%   Detailed explanation goes here

% tmp solution is that lines must be a matrix where each two consecutive
% rows are first x values in a odd column and the y values in the even
% column to the right, i.e.: 
% column 1 = line 1's x values
% column 2 = line 1's y values
% column 3 = line 2's x values
% column 4 = line 2's y values
% etc
% preferrably create a cell array with lines or some other more intuitive
% structure for lines input
    nLines = width(lines)/2;

    mse = 0;
    for i = 1:nLines
        line = lines(:, ((i*2)-1):(i*2));
        line = [line(~isnan(line(:,1)),1), line(~isnan(line(:,2)),2)];
        mse = mse + min(vecnorm(point - line, 2, 2))^2;
    end
    mse = mse/nLines;
    % dist1 = min(sqrt(sum((point - lines).^2, 2)));
    % dist2 = min(sqrt(sum((point - lines).^2, 2)));
    % dist3 = min(sqrt(sum((point - lines).^2, 2)));
    % dist4 = min(sqrt(sum((point - lines).^2, 2)));
    % mse = (dist1^2 + dist2^2  + dist3^2 + dist4^2)/4;%[dist1, dist2, dist3, dist4];inputArg2;
end

