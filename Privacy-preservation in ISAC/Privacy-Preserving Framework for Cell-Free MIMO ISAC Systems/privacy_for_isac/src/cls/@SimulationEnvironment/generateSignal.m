function [] = generateSignal(obj, nSymbols, modulationOrder)
%GENERATESIGNAL Summary of this function goes here
%   Detailed explanation goes here
    symbolIdx = randi([0, modulationOrder-1], [numel(obj.users), nSymbols]);
    symbols = qammod(symbolIdx, modulationOrder, 'UnitAveragePower', true);
    cNoise = sqrt(obj.noiseVar/2)*(randn(size(symbols)) + 1i*randn(size(symbols)));

    obj.cSignal = symbols + cNoise;

    %{
        num_waveforms = 1;
        codes = zeros(num_waveforms, nSymbols);
        for i = 1:num_waveforms
            codes(i,:) = sign(randn(1, nSymbols));
        end
        
        % Calculate cross-correlation matrix between codes
        %corr_matrix = codes * codes';
        
        % Uncorrelated waveforms using hadamard matrix (matrix with pairwise
        % orthogonal rows whos inner product is scaled N*I)
        % Generate uncorrelated waveforms using Hadamard matrix
        hadamard_matrix = hadamard(nSymbols);
        waveforms = zeros(num_waveforms, nSymbols);
        for i = 1:num_waveforms
            waveforms(i,:) = hadamard_matrix(i,:) * sqrt(nSymbols);
        end
        
        % Multiply each waveform with its corresponding code and normalize
        for i = 1:num_waveforms
            waveforms(i,:) = waveforms(i,:) .* codes(i,:) / norm(codes(i,:));
        end
    
        obj.sSignal = waveforms;
    %}

    % add noise time sqrt(2) since not complex symbols
    sNoise = sqrt(obj.noiseVar)*randn(1, nSymbols);
    obj.sSignal = (2*randi([0, 1], [1, nSymbols])-1) + sNoise;
    
    % TODO:
    % ('UnitAvaragePower', true) summan av alla Tx effekt = 1 kanske gångra
    % med 1/antal Tx?
end

