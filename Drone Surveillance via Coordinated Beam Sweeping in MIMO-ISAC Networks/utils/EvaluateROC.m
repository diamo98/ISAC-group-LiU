function [Pd_sorted_NC,Pfa_sorted_NC,Pd_sorted_WC,Pfa_sorted_WC,precoder_vec_pwr] = EvaluateROC(mod_order,MC,len_point,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03,pwr_bgt_perAP_db,Peta)

SSB_symbols =  qammod(randi([0,mod_order-1],MC,J,len_point),mod_order,UnitAveragePower=true); % [MC Q] can extend to [MC J Q]
Sens_symbols =  qammod(randi([0,mod_order-1],MC,len_point),mod_order,UnitAveragePower=true); % [MC Q]
pmin = max(max(rho_ssb_all)); 
%pwr_bgt_perAP_db = 20; 
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
T0_WC_sum_arr=zeros(MC,1);
T1_WC_sum_arr=zeros(MC,1);

for mc=1:MC % MC rounds at each point
    T0_WC_sum=0;
    T1_WC_sum=0;
    for q=1:len_point % Q iterate point
    %for mc=1:MC % MC rounds at each point
                if rcs_dis == 1 % Swerling 2 (flutuation every slot)
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
                    ( H01*(Sens_symbols(mc,q).*precoder_vec_pwr(:,1,q) + SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*conj(precoderSSB_GoB(:,q))) + ... %% direct link % both from SSB and Sensing (Sensing dominantate..)
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
                    ( H01*(Sens_symbols(mc,q).*precoder_vec_pwr(:,1,q) + SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*conj(precoderSSB_GoB(:,q))) + ... %% direct link % both from SSB and Sensing (Sensing dominantate..)
                    H02*(Sens_symbols(mc,q).*precoder_vec_pwr(:,2,q) + SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*conj(precoderSSB_GoB(:,q)) ) + ... 
                    H03*(Sens_symbols(mc,q).*precoder_vec_pwr(:,3,q) + SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*conj(precoderSSB_GoB(:,q)) ) ); 
                % combining V
                    y_beamformed_T0_WC = V(:,q)'*(y_received_T0);
                    y_beamformed_T1_WC = V(:,q)'*(y_received_T1);
                
                    y_beamformed_T0_NC = (y_received_T0);
                    y_beamformed_T1_NC = (y_received_T1);
                %end
                
                pwr_T0_WC = norm(y_beamformed_T0_WC)^2;
                pwr_T1_WC = norm(y_beamformed_T1_WC)^2; 
                pwr_T0_NC = norm(y_beamformed_T0_NC)^2;
                pwr_T1_NC = norm(y_beamformed_T1_NC)^2; 

                T0_WC_sum = T0_WC_sum+pwr_T0_WC;
                T1_WC_sum = T1_WC_sum+pwr_T1_WC; 

                T0_NC(q,mc) = pwr_T0_NC;
                T1_NC(q,mc) = pwr_T1_NC; 
    end % points
    T0_WC_sum_arr(mc)=T0_WC_sum;
    T1_WC_sum_arr(mc)=T0_WC_sum;
end % mc

% ---- ROC computation ----
% thresholds spanning both distributions
%Tmin_WC = min(min([T0_WC; T1_WC]));
%Tmax_WC = max(max([T0_WC; T1_WC]));
Tmin_WC=min([T0_WC_sum_arr; T1_WC_sum_arr]);
Tmax_WC=max([T0_WC_sum_arr; T1_WC_sum_arr]);
Tmin_NC = min(min([T0_NC; T1_NC]));
Tmax_NC = max(max([T0_NC; T1_NC]));
nTh  = 500;
Pd_NC = zeros(1,length(nTh));
Pfa_NC = zeros(1,length(nTh));
taus_NC = linspace(Tmin_NC, Tmax_NC, nTh);
Pd_WC = zeros(1,length(nTh));
Pfa_WC = zeros(1,length(nTh));
taus_WC = linspace(Tmin_WC, Tmax_WC, nTh);

for k = 1:nTh
    tau_WC = taus_WC(k);
    tau_NC = taus_NC(k);

    % Pd_arr_WC = sum(T1_WC > tau_WC);
    % Pd_WC(k)  = mean(Pd_arr_WC > 0);  % P_D = Pr(T > tau | H1) % one detect the whole sweep is detected
    % Pfa_arr_WC = sum(T0_WC > tau_WC);
    % Pfa_WC(k) = mean(Pfa_arr_WC > 0);  % P_FA = Pr(T > tau | H0)
    % 
    % Pd_arr_NC = sum(T1_NC > tau_NC); 
    % Pd_NC(k)  = mean(Pd_arr_NC > 0); % P_D = Pr(T > tau | H1) % one detect the whole sweep is detected
    % Pfa_arr_NC = sum(T0_NC > tau_NC); 
    % Pfa_NC(k) = mean(Pfa_arr_NC > 0); % P_FA = Pr(T > tau | H0)

     Pd_WC(k)  = (mean(T1_WC_sum_arr > tau_WC));  % P_D = Pr(T > tau | H1)
     Pfa_WC(k) = (mean(T0_WC_sum_arr > tau_WC));  % P_FA = Pr(T > tau | H0)
     
     Pd_NC(k)  = mean(mean(T1_NC > tau_NC));  % P_D = Pr(T > tau | H1)
     Pfa_NC(k) = mean(mean(T0_NC > tau_NC));  % P_FA = Pr(T > tau | H0)
end

[Pfa_sorted_NC, idx_NC] = sort(Pfa_NC, 'ascend');
Pd_sorted_NC = Pd_NC(idx_NC);

[Pfa_sorted_WC, idx_WC] = sort(Pfa_WC, 'ascend');
Pd_sorted_WC = Pd_WC(idx_WC);

end

