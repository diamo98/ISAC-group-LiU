function [] = optimizeConfiguration(obj, maxIts, threshold, maxPower, ueSinrConstraint)
%OPTIMIZECONFIGURATION Summary of this function goes here
%   Detailed explanation goes here
    tmpSimEnv = SimulationEnvironment(obj);
    %TODO: ugly to pick the first transmitAp to see how many antennas,  think of
    % something else maybe?
    nRAps = numel(obj.receiverAps);
    tmpReceiver = AccessPoint(obj.target.pos, obj.transmitAps(1).nAntennas);
    tmpSimEnv.transmitAps = [tmpSimEnv.transmitAps, tmpSimEnv.receiverAps];
    tmpSimEnv.receiverAps = tmpReceiver;
    tmpSimEnv.generateChannelProperties();

    % just so it doesnt break first iteration
    prevReceivers = obj.receiverAps;
    it = 0;
    % Extreme case when two APs are exatcly opposite e.g.(100, 50) and 
    % (50, 100) they would trigger while loop to stop?
    while ~all(ismember([prevReceivers.pos], [tmpSimEnv.receiverAps.pos]), 'all')
        it = it + 1;
        tmpSimEnv.computePrecoders(maxIts, threshold, maxPower, ueSinrConstraint);
    
        aggrSymbols = [tmpSimEnv.cSignal; tmpSimEnv.sSignal];
    
        % MIMO-bok sida 183 eq C.25
        %TODO: fixa index grejen snyggare med foreach
        mIs = [];
        for apIdx = 1:numel(tmpSimEnv.transmitAps)
            ueMis = zeros(1, numel(tmpSimEnv.users));
            for userIdx = 1:numel(tmpSimEnv.users)
                ap = tmpSimEnv.transmitAps(apIdx);
                x = ap.precoderMatrix * aggrSymbols;
                xs = ap.precoderMatrix(:, end)*tmpSimEnv.sSignal;
                userH = squeeze(tmpSimEnv.cH(userIdx, apIdx, :)).';
                %size(x) samma storlek
                %size(xs)
                %TODO: undersök snr konstanten? Inbakat i xs med precodern?
                ueMis(userIdx) = log2(1 + (1/tmpSimEnv.noiseVar)*userH*(xs*xs')*userH');
            end
            mIs = [mIs; real(ueMis)];
        end

        %For max sum of user MI per AP, ie row-wise sums sine its N_AP X N_UE
        %matrix
        mIs = sum(mIs, 2);
        sortedMIs = sort(mIs);
        sortedMIs = sortedMIs(end-nRAps+1:end);
        [~, I] = ismember(sortedMIs, mIs);
    
        prevReceivers = tmpSimEnv.receiverAps;
        tmpSimEnv.receiverAps = tmpSimEnv.transmitAps(I);
        tmpSimEnv.transmitAps(I) = [];
        tmpSimEnv.transmitAps = [tmpSimEnv.transmitAps, tmpSimEnv.receiverAps];
        tmpSimEnv.generateChannelProperties();
        if it > 50
            disp("Selection didnt converge");
            break;
        end
    end

    obj.receiverAps = tmpSimEnv.receiverAps;
    [~, I] = ismember(tmpSimEnv.receiverAps, tmpSimEnv.transmitAps);
    tmpSimEnv.transmitAps(I) = [];
    obj.transmitAps = tmpSimEnv.transmitAps;
    obj.generateChannelProperties();
end

