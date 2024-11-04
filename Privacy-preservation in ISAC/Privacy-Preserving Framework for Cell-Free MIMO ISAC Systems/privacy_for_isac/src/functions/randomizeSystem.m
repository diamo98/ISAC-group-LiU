function [usersLocs, apsLocs, targetLoc] = randomizeSystem(nUsers, nAps, gridSize)
% Generates locations for a system where the users are within the center of
% the accesspoints
    userRadius = gridSize * 0.60;
    userAngles = 2*pi*rand(1, nUsers);
    userAmp = userRadius*rand(1, nUsers);

    usersLocs = [floor(userAmp.*cos(userAngles))', floor(userAmp.*sin(userAngles))'];

    targetLoc = [floor(userRadius*rand()*cos(2*pi*rand)), floor(userRadius*rand()*sin(2*pi*rand))];

    apRadiusRange = [0.8*gridSize, gridSize];
    apAngles = zeros(1, nAps);
    apAmp = (apRadiusRange(2) - apRadiusRange(1))*rand(1, nAps) + apRadiusRange(1);
    apAngles(1) = 2*pi*rand();

    for i = 2:nAps
        % minimum spacing is a range percentage from equal range APs on the circle
        apAngles(i) = apAngles(i-1) + (2*pi/nAps)*((1.25-0.75)*rand() + 0.75);
    end

    apsLocs = [floor(apAmp.*cos(apAngles))', floor(apAmp.*sin(apAngles))'];
end

