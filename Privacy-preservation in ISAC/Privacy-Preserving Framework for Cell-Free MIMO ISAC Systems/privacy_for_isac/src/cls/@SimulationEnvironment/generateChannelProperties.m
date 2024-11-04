function [] = generateChannelProperties(obj, plExp)
%GENERATECHANNELPROPERTIES Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 1
        plExp = 2.8;
    end
    nUsers = numel(obj.users);
    nTAps = numel(obj.transmitAps);
    nRAps = numel(obj.receiverAps);
    nAntennas = obj.transmitAps.nAntennas;

    % communication
    cDists = zeros(nUsers, nTAps);
    for i = 1:nUsers
        for j = 1:nTAps
            cDists(i, j) = obj.users(i).distanceTo(obj.transmitAps(j));
        end
    end

    cBetas = (cDists.^-plExp).*ones(nUsers, nTAps);

    cH = sqrt(cBetas./2).*(randn(nUsers, nTAps, nAntennas) + 1j*randn(nUsers, nTAps, nAntennas));
    obj.cH = cH;
    obj.cBetas = cBetas;

    % sensing
    sDists = zeros(nTAps, nRAps);

    for i = 1:nTAps
        for j = 1:nRAps
            sDists(i,j) = obj.transmitAps(i).distanceTo(obj.target) + obj.receiverAps(j).distanceTo(obj.target);
        end
    end

    sBetas = (sDists.^-plExp).*ones(nTAps, nRAps);

    %sH = (1/sqrt(2))*(randn(nAntennas, nAntennas, nTAps) + 1j*randn(nAntennas, nAntennas, nTAps));
    %TODO: Blir typ en nRAps X nAntennas X nTAps X nAntennas, matris vilket
    %verkar hemskt hoppas det inte kommer behövas
    %sH = sqrt(cBetas./2).*(randn(nUsers, nTAps, nAntennas) + 1j*randn(nUsers, nTAps, nAntennas));
%     for i = 1:nTAps
%         sH(:,:,i) = sBetas(i)*sH(:,:,i);
%     end
% TODO:
% lös problemet

    obj.sH = nan;
    obj.sBetas = sBetas;
end

