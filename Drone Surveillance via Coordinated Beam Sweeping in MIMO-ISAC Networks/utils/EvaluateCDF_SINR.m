function [sorted_SCINR_db, cdf_vals,precoder_vec_pwr] = EvaluateCDF_SINR(mod_order,MC,len_point,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03)

% received signal y at APr [Mxq]
SSB_symbols =  qammod(randi([0,mod_order-1],MC,J,len_point),mod_order,UnitAveragePower=true); % [MC Q] can extend to [MC J Q]
Sens_symbols =  qammod(randi([0,mod_order-1],MC,len_point),mod_order,UnitAveragePower=true); % [MC Q]
pmin = min(min(rho_ssb_all)); 
pwr_bgt_perAP_db = 30; % 30 dBm
pwr_bgt_perAP = 10.^(pwr_bgt_perAP_db/10);
%weibull_lambda = 0.0400 ; weibull_k = 0.5; 
weibull_lambda = 0.0847;weibull_k = 0.833;
precoder_vec_pwr = zeros(M,J,len_point); % [M J Q]

for i_pwr = 1:length(pwr_bgt_perAP) % iterate over power level
    fprintf('Iterate Pmax: %d out of %d \n',i_pwr,length(pwr_bgt_perAP))
    Sens_power_budget_mat = pwr_bgt_perAP(i_pwr) - rho_ssb_all;
    Sens_power_budget_mat = (Sens_power_budget_mat>0).*Sens_power_budget_mat; % exclude negative power

    for k=1:len_point % allocate power accoding to each point 'q'
         precoder_vec_pwr(:,:,k) = precoder_vec(:,:,k)./vecnorm(precoder_vec(:,:,k), 2, 1) .* sqrt(Sens_power_budget_mat(:,k))';
    end
        SICNR_ite = zeros(len_point,MC);
    for q=1:len_point % Q iterate point
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
            
            y_received_numer = Sens_symbols(mc,q).*( ... 
                rcs1.*H_all_fullloss(:,1:M,q)*precoder_vec_pwr(:,1,q) + ... 
                rcs2.*H_all_fullloss(:,M+1:2*M,q)*precoder_vec_pwr(:,2,q) + ... 
                rcs3.*H_all_fullloss(:,2*M+1:3*M,q)*precoder_vec_pwr(:,3,q));

            y_received_nu_pwr = abs(V(:,q)'*y_received_numer)^2;
            
            y_received_denome = sqrt(beta_clutter/2).*( ... 
                SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*(randn(M,M)+1i*randn(M,M)) +  ...
                SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*(randn(M,M)+1i*randn(M,M)) + ...
                SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*(randn(M,M)+1i*randn(M,M)) ) ...
                *conj(precoderSSB_GoB(:,q))   + sqrt(noise_var/2).*(randn(M,1)+1i*randn(M,1)) + ... 
              ( H01*(Sens_symbols(mc,q).*precoder_vec_pwr(:,1,q) + SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*conj(precoderSSB_GoB(:,q))) + ... 
                H02*(Sens_symbols(mc,q).*precoder_vec_pwr(:,2,q) + SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*conj(precoderSSB_GoB(:,q)) ) + ... 
                H03*(Sens_symbols(mc,q).*precoder_vec_pwr(:,3,q) + SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*conj(precoderSSB_GoB(:,q)) ) ); %% direct link
            
            y_received_de_comb_pwr = abs(V(:,q)'*y_received_denome)^2;
            
            SICNR_ite(q,mc) = y_received_nu_pwr/y_received_de_comb_pwr;
        end

    end % iterate i (point q)
end % for pwr arr

SICNR_array = reshape(SICNR_ite,1,[]);

sorted_SCINR = sort(SICNR_array);
sorted_SCINR_db = 10*log10(sorted_SCINR);

% Compute empirical CDF
n = length(sorted_SCINR);
cdf_vals = (1:n) / n;



end