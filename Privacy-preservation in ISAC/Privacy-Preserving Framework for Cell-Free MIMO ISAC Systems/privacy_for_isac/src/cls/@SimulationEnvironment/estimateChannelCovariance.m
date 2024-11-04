function [] = estimateChannelCovariance(obj, maxIts, threshold, channelEstVar)
%ESTIMATECHANNELCOVARIANCE Summary of this function goes here
%   Detailed explanation goes here

    % TODO: 
    % 1. Ta in mottagna signalen som argument?
    % 2. Deep expectation maximization artikel är H en reell version av den
    % komplexa kanalmatrisen som ser ut på ett speciellt sätt men som inte
    % har brytts om här?
    % 3.
    hEstVar = channelEstVar;%1e-10;%1e-2; % gör om till argument
    dOF = 10; % degrees of freedom, gör om till argument.
    % flytta till User? :/

    M = obj.transmitAps.nAntennas;
    % TODO: 
    % ha en user som är adversary istället för bara sista indexet
    %adversaryUserIdx = find(obj.users == obj.adversary);
    totR = zeros(M, M, numel(obj.transmitAps));
    nSymbols = width(obj.sSignal);


    for curTApIdx = 1:numel(obj.transmitAps)
        W = obj.transmitAps(curTApIdx).precoderMatrix;
        W(:,obj.adversaryIdx) = [];
        R_tmp = cell(1, nSymbols);
        % tmp populate cov cellarray
        for i = 1:nSymbols
            R_tmp{i} = zeros(M);
        end
        for j = 1:(obj.users(obj.adversaryIdx).nAntennas * obj.nObservations)
            obj.generateChannelProperties();
            for i = 1:nSymbols
    
                % TODO: What is happening with end-1??????? Världens
                % fiasko, men kanalerna görs separat för kommunikation och
                % sensing så -1 måste varit en fulfix för att
                % kommunikationskanalkeofficienterna skulle funka för både
                % (nästan) alla kommunikationssignalr same sensing
                % signalen...
                % Eller för att man antas ta bort sin egen signal som
                % adversary ur precodern????? För funkar i optimize?
                aggrSymbols = [obj.cSignal(1:end-1,i); obj.sSignal(i)];
                adversaryH = squeeze(obj.cH(obj.adversaryIdx, curTApIdx, :)).';
                %adversaryH = adversaryH
                y = adversaryH*W*aggrSymbols;
    
                hEst = adversaryH + sqrt((hEstVar)/2)*(randn(size(adversaryH)) + 1i*randn(size(adversaryH)));
                
                % initial x estimate (zero forcing)
                % TODO: kika på varningens orsak
                warning("off", "MATLAB:nearlySingularMatrix");
                xEst = (hEst'*hEst)\hEst'*y;
                warning("on", "MATLAB:nearlySingularMatrix");
    
                % initial E(H)
                eHEst = hEst;
                %eHEst = adversaryH + sqrt((1e-2)/2)*(randn(size(adversaryH)) + 1i*randn(size(adversaryH)));
    
                % initial omega
                omega = randn(size(M)) + 1i*randn(size(M));
                for k = 1:maxIts
                    % calculate q(u) and E(u)
                    xEstPrev = xEst;
                    % TODO: dOF borde också vara en estimerad parameter?
                    alpha = (dOF + 2*obj.users(obj.adversaryIdx).nAntennas)/2; % 50% med 2*M_r reciver antennas
                    C = 1/obj.noiseVar*((y-eHEst*xEst)'*(y-eHEst*xEst) + xEst'*(M*omega)*xEst); % missing term?
                    beta = (dOF + C)/2;
                    qApprx = gamrnd(alpha, beta);
                    uMean = alpha/beta;
    
                    % calculate omega and E(q(H))
                    omega = inv((1/hEstVar)*eye(M) + (uMean/obj.noiseVar)*(xEst*xEst'));
                    eHEst = (omega*((1/hEstVar)*hEst' + (uMean/obj.noiseVar)*(xEst*y')))';
                    V = chol(real(M*omega));
                    %xEst = min(norm(y-eHEst*xEst)^2, norm(V*xEst)^2);  eq 29?
    
                    % find x as minimum to eq 29
                    %opFun = @(x) norm(y-eHEst*x)^2 + norm(V*x)^2;
                    %xEst = fminsearch(opFun, xEst, optimset('Display','off'));
                    Htot = [eHEst;V];
                    yTot = [y; zeros(M, 1)];
                    xEst = (Htot'*Htot)\Htot'*yTot;
                    %hEst = eHEst; % hEst updateras aldrig i eHEst
                    %uttrycket

                    if norm(xEst - xEstPrev) <= threshold
                        break;
                    end
                end
                R_tmp{i} = R_tmp{i} + xEst*xEst';
            end
        end
        R = R_tmp;
        obj.estimatedR{curTApIdx} = R;
    end
    %obj.estimatedR{1}{2}
    %obj.estimatedR = totR;
end

