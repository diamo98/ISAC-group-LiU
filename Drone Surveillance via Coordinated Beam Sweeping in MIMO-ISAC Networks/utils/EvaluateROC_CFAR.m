function [Pd_sorted_NC,Pfa_sorted_NC,Pd_sorted_WC,Pfa_sorted_WC,precoder_vec_pwr, ...
          Pd_sorted_CFAR_NC,Pfa_sorted_CFAR_NC,Pd_sorted_CFAR_WC,Pfa_sorted_CFAR_WC] = ...
    EvaluateROC_CFAR(mod_order,MC,len_point,rho_ssb_all,rcs_var_sim, ...
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03,pwr_bgt_perAP_db,Peta)

% -----------------------------
% (unchanged initial part)
% -----------------------------
SSB_symbols =  qammod(randi([0,mod_order-1],MC,J,len_point),mod_order,UnitAveragePower=true); % [MC Q] can extend to [MC J Q]
Sens_symbols =  qammod(randi([0,mod_order-1],MC,len_point),mod_order,UnitAveragePower=true); % [MC Q]
pmin = max(max(rho_ssb_all)); 
pwr_bgt_perAP = 10.^(pwr_bgt_perAP_db/10);
Sens_power_budget_mat = pwr_bgt_perAP - rho_ssb_all; %

precoder_vec_pwr = zeros(M,J,len_point);
for k=1:len_point % allocate power accoding to each point 'q'
    precoder_vec_pwr(:,:,k) = precoder_vec(:,:,k)./vecnorm(precoder_vec(:,:,k), 2, 1) .* sqrt(Sens_power_budget_mat(:,k))';
end

T1_NC = zeros(len_point,MC);    % test statistics under H1
T0_NC = zeros(len_point,MC);    % test statistics under H0
T1_WC = zeros(len_point,MC);    % test statistics under H1
T0_WC = zeros(len_point,MC);    % test statistics under H0

for q=1:len_point % Q iterate point
    for mc=1:MC % MC rounds at each point
        if rcs_dis == 1 % Swerling 2 (fluctuate every slot)
            rcs1=sqrt(rcs_var_sim/2).*(randn()+1i*randn());
            rcs2=sqrt(rcs_var_sim/2).*(randn()+1i*randn());
            rcs3=sqrt(rcs_var_sim/2).*(randn()+1i*randn()); 
        elseif rcs_dis == 2 % Weibull distribution
            rcs1= sqrt(wblrnd(weibull_lambda, weibull_k, [M, 1]));
            rcs2= sqrt(wblrnd(weibull_lambda, weibull_k, [M, 1]));
            rcs3= sqrt(wblrnd(weibull_lambda, weibull_k, [M, 1]));
        end
        clutter_reli1 = (randn(M,M)+1i*randn(M,M));
        clutter_reli2 = (randn(M,M)+1i*randn(M,M));
        clutter_reli3 = (randn(M,M)+1i*randn(M,M));
        AWGN = sqrt(noise_var/2).*(randn(M,1)+1i*randn(M,1)); % AWGN

        y_received_T0 = 0 + ...
            sqrt(beta_clutter/2).*( ...  % clutter from SSB
            SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*clutter_reli1 +  ...
            SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*clutter_reli2 + ...
            SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*clutter_reli3 )*conj(precoderSSB_GoB(:,q)) ...
            + AWGN + ... % AWGN
            ( H01*(Sens_symbols(mc,q).*precoder_vec_pwr(:,1,q) + SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*conj(precoderSSB_GoB(:,q))) + ... %% direct link
            H02*(Sens_symbols(mc,q).*precoder_vec_pwr(:,2,q) + SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*conj(precoderSSB_GoB(:,q)) ) + ... 
            H03*(Sens_symbols(mc,q).*precoder_vec_pwr(:,3,q) + SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*conj(precoderSSB_GoB(:,q)) ) );

        eta_n = rand(1)<Peta;
        y_received_T1 = eta_n.*Sens_symbols(mc,q).*( ... % sensing echoes [M x 1]
            rcs1.*H_all_fullloss(:,1:M,q)*precoder_vec_pwr(:,1,q) + ... 
            rcs2.*H_all_fullloss(:,M+1:2*M,q)*precoder_vec_pwr(:,2,q) + ... 
            rcs3.*H_all_fullloss(:,2*M+1:3*M,q)*precoder_vec_pwr(:,3,q)) + ...
            sqrt(beta_clutter/2).*( ...  % clutter from SSB
            SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*clutter_reli1 +  ...
            SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*clutter_reli2 + ...
            SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*clutter_reli3 )*conj(precoderSSB_GoB(:,q)) ...
            + AWGN + ... % AWGN
            ( H01*(Sens_symbols(mc,q).*precoder_vec_pwr(:,1,q) + SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*conj(precoderSSB_GoB(:,q))) + ... %% direct link
            H02*(Sens_symbols(mc,q).*precoder_vec_pwr(:,2,q) + SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*conj(precoderSSB_GoB(:,q)) ) + ... 
            H03*(Sens_symbols(mc,q).*precoder_vec_pwr(:,3,q) + SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*conj(precoderSSB_GoB(:,q)) ) ); 

        % combining V
        y_beamformed_T0_WC = V(:,q)'*(y_received_T0);
        y_beamformed_T1_WC = V(:,q)'*(y_received_T1);

        y_beamformed_T0_NC = (y_received_T0);
        y_beamformed_T1_NC = (y_received_T1);

        pwr_T0_WC = norm(y_beamformed_T0_WC)^2;
        pwr_T1_WC = norm(y_beamformed_T1_WC)^2; 
        pwr_T0_NC = norm(y_beamformed_T0_NC)^2;
        pwr_T1_NC = norm(y_beamformed_T1_NC)^2; 

        T0_WC(q,mc) = pwr_T0_WC;
        T1_WC(q,mc) = pwr_T1_WC; 
        T0_NC(q,mc) = pwr_T0_NC;
        T1_NC(q,mc) = pwr_T1_NC; 
    end
end

% ---- ROC computation for energy detector (unchanged) ----
Tmin_WC = min(min([T0_WC; T1_WC]));
Tmax_WC = max(max([T0_WC; T1_WC]));
Tmin_NC = min(min([T0_NC; T1_NC]));
Tmax_NC = max(max([T0_NC; T1_NC]));
nTh  = 500;
Pd_NC = zeros(1,nTh);
Pfa_NC = zeros(1,nTh);
taus_NC = linspace(Tmin_NC, Tmax_NC, nTh);
Pd_WC = zeros(1,nTh);
Pfa_WC = zeros(1,nTh);
taus_WC = linspace(Tmin_WC, Tmax_WC, nTh);

for k = 1:nTh
    tau_WC = taus_WC(k);
    tau_NC = taus_NC(k);

    Pd_WC(k)  = mean(mean(T1_WC > tau_WC));  % P_D = Pr(T > tau | H1)
    Pfa_WC(k) = mean(mean(T0_WC > tau_WC));  % P_FA = Pr(T > tau | H0)
    Pd_NC(k)  = mean(mean(T1_NC > tau_NC));  % P_D = Pr(T > tau | H1)
    Pfa_NC(k) = mean(mean(T0_NC > tau_NC));  % P_FA = Pr(T > tau | H0)
end

[Pfa_sorted_NC, idx_NC] = sort(Pfa_NC, 'ascend');
Pd_sorted_NC = Pd_NC(idx_NC);

[Pfa_sorted_WC, idx_WC] = sort(Pfa_WC, 'ascend');
Pd_sorted_WC = Pd_WC(idx_WC);

% -----------------------------
% CFAR (CA-CFAR) ROC computation
% -----------------------------
% CFAR parameters (tune these as desired)
num_train = 8;    % must be even (total training cells)
num_guard = 2;    % total guard cells (usually even)
if mod(num_train,2)~=0
    error('num_train must be even for symmetric windowing');
end
half_train = num_train/2;
half_guard = num_guard/2;

% vector of desired Pfa values to sweep for ROC
Pfa_vec = logspace(-5,-1,50);  % change as needed
nPfa = length(Pfa_vec);

Pd_cfar_NC = zeros(1,nPfa);
Pfa_cfar_NC = zeros(1,nPfa);
Pd_cfar_WC = zeros(1,nPfa);
Pfa_cfar_WC = zeros(1,nPfa);

% we'll only test CFAR for indices where full window fits:
start_idx = 1 + half_train + half_guard;
end_idx   = len_point - (half_train + half_guard);
valid_idx = start_idx:end_idx;
num_valid = length(valid_idx);

% If there are no valid indices, CFAR cannot be applied:
if num_valid < 1
    warning('len_point too small for the chosen CFAR window. CFAR ROC will be empty.');
    Pd_sorted_CFAR_NC = [];
    Pfa_sorted_CFAR_NC = [];
    Pd_sorted_CFAR_WC = [];
    Pfa_sorted_CFAR_WC = [];
    return;
end

% for each desired Pfa compute CA-CFAR detection across q (range cells) and MC
for p_idx = 1:nPfa
    Pfa_des = Pfa_vec(p_idx);
    % scale factor alpha for CA-CFAR (exponential noise): alpha = N*(Pfa^{-1/N} - 1)
    N = num_train;
    alpha = N*(Pfa_des^(-1/N) - 1);

    % detection decisions containers
    D0_NC = false(num_valid, MC); % H0 decisions
    D1_NC = false(num_valid, MC); % H1 decisions
    D0_WC = false(num_valid, MC);
    D1_WC = false(num_valid, MC);

    for mc = 1:MC
        % build vectors along q for this Monte Carlo realization
        vec_T0_NC = T0_NC(:,mc); % len_point x 1
        vec_T1_NC = T1_NC(:,mc);
        vec_T0_WC = T0_WC(:,mc);
        vec_T1_WC = T1_WC(:,mc);

        for ii = 1:num_valid
            q = valid_idx(ii);

            % training indices: left and right excluding guard cells and CUT
            left_train = (q - half_guard - half_train):(q - half_guard - 1);
            right_train = (q + half_guard + 1):(q + half_guard + half_train);

            training_idx = [left_train, right_train];

            % noise estimate = average power in training cells (cell-averaging)
            noise_est_NC = mean(vec_T0_NC(training_idx)); % under H0 this is appropriate
            noise_est_WC = mean(vec_T0_WC(training_idx));

            % threshold
            thr_NC = alpha * noise_est_NC;
            thr_WC = alpha * noise_est_WC;

            % compare with CUT for both H0 and H1 (this yields detection decision)
            D0_NC(ii,mc) = vec_T0_NC(q) > thr_NC;
            D1_NC(ii,mc) = vec_T1_NC(q) > thr_NC; % use same thr_NC estimated from local H0 estimate
            D0_WC(ii,mc) = vec_T0_WC(q) > thr_WC;
            D1_WC(ii,mc) = vec_T1_WC(q) > thr_WC;
        end
    end

    % compute Pfa & Pd as averages across all valid CUT positions and MC trials
    Pfa_cfar_NC(p_idx) = mean(D0_NC(:));
    Pd_cfar_NC(p_idx)  = mean(D1_NC(:));
    Pfa_cfar_WC(p_idx) = mean(D0_WC(:));
    Pd_cfar_WC(p_idx)  = mean(D1_WC(:));
end

% sort outputs for ROC plotting
[Pfa_sorted_CFAR_NC, idx_cfarnc] = sort(Pfa_cfar_NC, 'ascend');
Pd_sorted_CFAR_NC = Pd_cfar_NC(idx_cfarnc);

[Pfa_sorted_CFAR_WC, idx_cfarwc] = sort(Pfa_cfar_WC, 'ascend');
Pd_sorted_CFAR_WC = Pd_cfar_WC(idx_cfarwc);

end
