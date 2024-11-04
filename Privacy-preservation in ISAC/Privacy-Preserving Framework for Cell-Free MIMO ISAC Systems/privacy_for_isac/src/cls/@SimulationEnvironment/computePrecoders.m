function [] = computePrecoders(obj, maxIts, threshold, maxPower, ueSinrConstraint)
%COMPUTEPRECODERS Summary of this function goes here
%   Detailed explanation goes here

% techincally these precoder computations assume perfect CSI since the real
% channel is used for the ueSinr constraint?


    % max power per ap
    maxP = maxPower;%10^(50/10)/1000; % from isabellas, 50dBm

    %threshold for the communication
    sinrUeMin = ueSinrConstraint;%10^(3/10); % from isabellas, 3dB

    nUsers = numel(obj.users);
    M = obj.transmitAps.nAntennas;
    nTAps = numel(obj.transmitAps);
    nRAps = numel(obj.receiverAps);

    % total transmission signal
    signal = [obj.sSignal;
              obj.cSignal];

    % intial precoder (without receiver AP)
    % sensing precoder is last column
    unsolvableCnt = 0;
    W = randn(M*nTAps,nUsers+1)+1i*randn(M*nTAps,nUsers+1);
    %W = -20 + (20 + 20)*rand(M*nTAps,nUsers+1);

    % construct concatenated A matrix
    ATot = zeros(M*(nTAps));
    %TODO: fix names and make this multiple receiver aps patch less bad...
    ATotTot = zeros([size(ATot), nRAps]);
    % check every pair of APs like described in the thesis
 
 

    for k = 1:nRAps
        for i = 1:nTAps
            for j = 1:nTAps
                % steering vectors: AP1, AP2, and receiver
                % plus pi since supposed to be from target to APs according to paper?
                a1 = obj.transmitAps(i).calcSteeringVector(obj.target, "from");
                a2 = obj.transmitAps(j).calcSteeringVector(obj.target, "from");
                ar = obj.receiverAps(k).calcSteeringVector(obj.target, "from");
    
                A = sqrt(obj.sBetas(i, k)*obj.sBetas(j, k))*...
                    conj(a1)*ar'*(obj.target.rcs(i)*conj(obj.target.rcs(j)))*ar*a2.';
                ATot(((i-1)*M+1):(i*M),((j-1)*M+1):(j*M)) = A;
            end
        end
        ATotTot(:,:,k) = ATot;
    end

    %cvx_optval
    for i = 1:maxIts
        prevW = W;
        cvx_clear
        cvx_begin quiet

            variable W((M*(nTAps)), (nUsers+1)) complex
            variable tau_n nonnegative
            variable tau_d nonnegative

            % objective function
            % sinrSTot = 0;
            % for j = 1:width(obj.cSignal)
            %     sinrSTot = sinrSTot + real(signal(:,j)'*...
            %         (2*prevW'*Atot*(W-prevW) + prevW'*Atot*prevW)...
            %         *signal(:,j));
            % end
            
            % sinrSTot = sinrSTot / (width(signal)*M*obj.noiseVar);
            sinrSTot = 0;
            for j = 1:nRAps                
                sinrSTot = sinrSTot + trace(real(signal'*(2*prevW'*ATotTot(:,:,j)*(W-prevW) + prevW'*ATotTot(:,:,j)*prevW)*signal)...
                    /(width(signal)*nRAps*M*obj.noiseVar));
            end

            maximize sinrSTot
            subject to
                % User SINR constraint
                for k = 1:nUsers
                    hUser = reshape(obj.cH(k,:,:),1,[])';

                    otherUsers = 1:nUsers;
                    otherUsers(k) = [];

                    wUser = W(:,k);
                    wOtherUsers = W(:,otherUsers);
                    wSensing = W(:,end);
                    wUserPrev = prevW(:,k);

                    num = real(wUserPrev'*(hUser*hUser')*wUserPrev+...
                        2*wUserPrev'*(hUser*hUser')*(wUser-wUserPrev));
                    if isempty(otherUsers)
                        denom = pow_abs(hUser'*wSensing, 2) +...
                        obj.noiseVar;
                    else
                        denom = pow_abs(sum(hUser'*wOtherUsers), 2) +...
                        pow_abs(hUser'*wSensing, 2) +...
                        obj.noiseVar;
                    end

                    num >= tau_n;
                    denom <= tau_d;
                    tau_n >= tau_d*sinrUeMin;
                end

                % Power constraint
                for k = 1:nTAps
                    apW = W((((k-1)*M)+1):(k*M), :);
                    sum(sum(pow_abs(apW, 2), 2)) <= maxP;
                    % for l = 1:M
                    %     sum(pow_abs(apW(l,:), 2)) <= maxP;
                    % end
                end
            cvx_end

        if isnan(W)
            % sätt ett tak för o-konvergens iterationer och slumpa om
            % precodermatrisen
            unsolvableCnt = unsolvableCnt + 1;
            if unsolvableCnt == maxIts
                disp("precoder didnt converge");
            end
            % W = (1/sqrt(unsolvableCnt))*randn(M*nTAps,nUsers+1)+1i*randn(M*nTAps,nUsers+1);
            W = randn(M*nTAps,nUsers+1)+1i*randn(M*nTAps,nUsers+1);
        end

        if norm(W - prevW, 'fro')/norm(W, 'fro') < threshold
            break;
        end
    end

    for i = 1:nTAps
        apW = W((((i-1)*M)+1):(i*M), :);
        obj.transmitAps(i).precoderMatrix = apW;
    end
    obj.precoders = W;
end

