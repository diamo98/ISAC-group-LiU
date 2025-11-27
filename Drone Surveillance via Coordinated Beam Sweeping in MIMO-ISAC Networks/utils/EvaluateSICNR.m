function [SICNR_sim_mean,SICNR_exp_mean,pwr_bgt_perAP_db,SICNR_sim_mean_uncoordi] = EvaluateSICNR(mod_order,MC,len_point,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03,uncoordi_prec_all)

% received signal y at APr [Mxq]
SSB_symbols =  qammod(randi([0,mod_order-1],MC,J,len_point),mod_order,UnitAveragePower=true); % [MC Q] can extend to [MC J Q]
Sens_symbols =  qammod(randi([0,mod_order-1],MC,len_point),mod_order,UnitAveragePower=true); % [MC Q]
pmin = min(min(rho_ssb_all)); pmin_dbm = 10*log10(pmin);
pwr_bgt_perAP_db = [pmin_dbm:0.1:ceil(pmin_dbm)-0.1, ceil(pmin_dbm):0.5:30]; % power budget % MODIFY ME 1 to 0.5
pwr_bgt_perAP = 10.^(pwr_bgt_perAP_db/10);
len_pwr_budgt = length(pwr_bgt_perAP);
%(min(pwr_bgt_perAP) > min(min(rho_ssb_all)))  % sanity check
SICNR_sim_mat = zeros(len_point,len_pwr_budgt); % [Q pwr_len] precoder_vec_pwr = (precoder_vec); %
SICNR_sim_mat_uncoordi = zeros(len_point,len_pwr_budgt);
SICNR_exp_mat =  zeros(len_point,len_pwr_budgt); % [Q pwr_len] 
%weibull_lambda = 0.0400 ; weibull_k = 0.5; 
weibull_lambda = 0.0847;weibull_k = 0.833;
%precoder_vec_pwr = zeros(M,J,len_point); % [M J Q]

for i_pwr = 1:len_pwr_budgt % iterate over power level
    fprintf('Iterate Pmax: %d out of %d \n',i_pwr,len_pwr_budgt)
    Sens_power_budget_mat = pwr_bgt_perAP(i_pwr) - rho_ssb_all;
    Sens_power_budget_mat = (Sens_power_budget_mat>0).*Sens_power_budget_mat; % exclude negative power
    %for k=1:len_point % allocate power accoding to each point 'q'
    %     precoder_vec_pwr(:,:,k) = precoder_vec(:,:,k)./vecnorm(precoder_vec(:,:,k), 2, 1) .* Sens_power_budget_mat(:,k)';
    %end
    for q=1:len_point % Q iterate point
        y_received_nu_comb = 0; %zeros(M,len_mat,MC);
        y_received_denom_pwr_comb = 0; % zeros(M,len_mat,MC);
        y_received_nu_noncoordi_comb = 0; %zeros(M,len_mat,MC);
        y_received_denom_pwr_noncoordi_comb = 0; % zeros(M,len_mat,MC);

        precoder_vec_AP1=precoder_vec(:,1,q)./norm(precoder_vec(:,1,q)).*sqrt(Sens_power_budget_mat(1,q));
        precoder_vec_AP2=precoder_vec(:,2,q)./norm(precoder_vec(:,2,q)).*sqrt(Sens_power_budget_mat(2,q));
        precoder_vec_AP3=precoder_vec(:,3,q)./norm(precoder_vec(:,3,q)).*sqrt(Sens_power_budget_mat(3,q));
        precoder_vec_AP1_noncoordi=uncoordi_prec_all(:,1,q)./norm(uncoordi_prec_all(:,1,q)).*sqrt(Sens_power_budget_mat(1,q)); % noncoordinated precoder
        precoder_vec_AP2_noncoordi=uncoordi_prec_all(:,2,q)./norm(uncoordi_prec_all(:,2,q)).*sqrt(Sens_power_budget_mat(2,q));
        precoder_vec_AP3_noncoordi=uncoordi_prec_all(:,3,q)./norm(uncoordi_prec_all(:,3,q)).*sqrt(Sens_power_budget_mat(3,q));
        
        for mc=1:MC
            if rcs_dis == 1 % Swerling 2 (flutuation every slot)
                rcs1=sqrt(rcs_var_sim/2).*(randn()+1i*randn());
                rcs2=sqrt(rcs_var_sim/2).*(randn()+1i*randn());
                rcs3=sqrt(rcs_var_sim/2).*(randn()+1i*randn()); 
            elseif rcs_dis == 2 % Weibull distribution
                rcs1= sqrt(wblrnd(weibull_lambda, weibull_k, [M, 1]));
                rcs2= sqrt(wblrnd(weibull_lambda, weibull_k, [M, 1]));
                rcs3= sqrt(wblrnd(weibull_lambda, weibull_k, [M, 1]));
            end
            % proposed precoder 
            y_received_numer = Sens_symbols(mc,q).*( ...  % 1st term in (6)
                rcs1.*H_all_fullloss(:,1:M,q)*precoder_vec_AP1 + ... 
                rcs2.*H_all_fullloss(:,M+1:2*M,q)*precoder_vec_AP2 + ... 
                rcs3.*H_all_fullloss(:,2*M+1:3*M,q)*precoder_vec_AP3 );
            y_received_nu_pwr = abs(V(:,q)'*y_received_numer)^2; %combiner
            %y_received_numer_n = y_received_numer./norm(y_received_numer);
            %y_received_nu_pwr = abs(y_received_numer_n'*y_received_numer)^2; % optimial detector

            y_received_nu_comb = y_received_nu_comb + y_received_nu_pwr;
            y_received_denome = sqrt(beta_clutter/2).*( ...  % 2nd term in (6)
                SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*(randn(M,M)+1i*randn(M,M)) +  ...
                SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*(randn(M,M)+1i*randn(M,M)) + ...
                SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*(randn(M,M)+1i*randn(M,M)) )*conj(precoderSSB_GoB(:,q))  + ...
                  sqrt(noise_var/2).*(randn(M,1)+1i*randn(M,1)) + ...  % 4th term in (6)
               (H01*(Sens_symbols(mc,q).*precoder_vec_AP1 + SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*conj(precoderSSB_GoB(:,q)) ) + ... % 3rd term in (6)
                H02*(Sens_symbols(mc,q).*precoder_vec_AP2 + SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*conj(precoderSSB_GoB(:,q)) ) + ... 
                H03*(Sens_symbols(mc,q).*precoder_vec_AP3 + SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*conj(precoderSSB_GoB(:,q)) ) ); %% direct link
            y_received_de_pwr = abs(V(:,q)'*y_received_denome)^2;
            %y_received_de_pwr = abs(y_received_numer_n'*y_received_denome)^2;
            y_received_denom_pwr_comb = y_received_denom_pwr_comb + y_received_de_pwr;

            % non coordinated beamformer
            y_received_numer_noncoordi = Sens_symbols(mc,q).*( ...  % 1st term in (6)
                rcs1.*H_all_fullloss(:,1:M,q)*precoder_vec_AP1_noncoordi + ... 
                rcs2.*H_all_fullloss(:,M+1:2*M,q)*precoder_vec_AP2_noncoordi + ... 
                rcs3.*H_all_fullloss(:,2*M+1:3*M,q)*precoder_vec_AP3_noncoordi );
            y_received_nu_pwr_noncoordi = abs(V(:,q)'*y_received_numer_noncoordi)^2;
            y_received_nu_noncoordi_comb = y_received_nu_noncoordi_comb + y_received_nu_pwr_noncoordi;
            y_received_denome_noncoordi = sqrt(beta_clutter/2).*( ...  % 2nd term in (6)
                SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*(randn(M,M)+1i*randn(M,M)) +  ...
                SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*(randn(M,M)+1i*randn(M,M)) + ...
                SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*(randn(M,M)+1i*randn(M,M)) )*conj(precoderSSB_GoB(:,q))  + ...
                  sqrt(noise_var/2).*(randn(M,1)+1i*randn(M,1)) + ...  % 4th term in (6)
               (H01*(Sens_symbols(mc,q).*precoder_vec_AP1_noncoordi + SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*conj(precoderSSB_GoB(:,q)) ) + ... % 3rd term in (6)
                H02*(Sens_symbols(mc,q).*precoder_vec_AP2_noncoordi + SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*conj(precoderSSB_GoB(:,q)) ) + ... 
                H03*(Sens_symbols(mc,q).*precoder_vec_AP3_noncoordi + SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*conj(precoderSSB_GoB(:,q)) ) ); %% direct link
            y_received_de_pwr_noncoordi = abs(V(:,q)'*y_received_denome_noncoordi)^2;
            y_received_denom_pwr_noncoordi_comb = y_received_denom_pwr_noncoordi_comb + y_received_de_pwr_noncoordi;
            
        end
        
        SICNR_sim_mat(q,i_pwr) = y_received_nu_comb/y_received_denom_pwr_comb;
        SICNR_sim_mat_uncoordi(q,i_pwr) = y_received_nu_noncoordi_comb/y_received_denom_pwr_noncoordi_comb;

        y_received_nu_exp = rcs_var_sim*( abs(V(:,q)'*H_all_fullloss(:,1:M,q)*precoder_vec_AP1)^2 + abs(V(:,q)'*H_all_fullloss(:,M+1:2*M,q)*precoder_vec_AP2)^2 +...
                            abs(V(:,q)'*H_all_fullloss(:,2*M+1:3*M,q)*precoder_vec_AP3)^2 ); 
        y_received_denom_exp = sum(rho_ssb_all(:,q))*beta_clutter + noise_var;
        SICNR_exp_mat(q,i_pwr) = y_received_nu_exp/y_received_denom_exp;
    
    end % iterate i (point q)
end % for pwr arr

SICNR_sim_mean = mean(SICNR_sim_mat,1);
SICNR_sim_mean_uncoordi = mean(SICNR_sim_mat_uncoordi,1);
SICNR_exp_mean = mean(SICNR_exp_mat,1);


end