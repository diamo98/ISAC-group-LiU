%% Initial Setup %%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%% setting: %%%%%%%%%%%%%
addpath('./utils');
rcs_var_db = -10; gamma_ue_db = 3 ; %%
z=10; Pmax = 10^(30/10); %%%%% MODIFY PARAMETERS here!
MV=12; MH=12; M=MV*MH; J=3; plot_flag = 1; case_interf = 1; % 0: no freq reuse, freq reuse
% %%%%% r-radius z-level, d-spacing, r_outer-hexagonal
r = 250; d=2; rcs_var = 10^(rcs_var_db/10); noise_var = 10^(-60/10);   
gamma_ue =  10^(gamma_ue_db/10); fc=15e9; lambda=3e8/fc; r_outer = r/cosd(30); % hexagonal: r is inner radius and r_outer is outer radius
% %%%%% define a space volumn [width x height x depth]
vol_width = 6 ; vol_height = 2; vol_depth = 2;
point_row=ceil(vol_height/d +1); point_col=ceil(vol_width/d +1); point_deep=ceil(vol_depth/d +1);
tot_points = point_row*point_col*point_deep; mod_order = 4;

% define cartisian of points and APs 
[AP_positions,point_positions] = setup_cart(r,z,point_row,point_col,point_deep,d,plot_flag);
% define polar between APs_tx and points % azi_direct direct link
[azi_deg_tx2points, ele_deg_tx2points, radius_tx2points, azi_direct_dep, azi_direct_ari] = polar_APs_points(AP_positions,point_positions); 
% define polar between points and APs_rx
[azi_deg_points2rx, ele_deg_points2rx, radius_points2rx] = polar_Points_APrx(AP_positions,point_positions);

% W= [W1,W2,...,W_R] SSB beams and 
[azimu_GoB_SSB_deg, polar_GoB_SSB_deg, BeamList, nr_tot_SSBbeams, precoderSSB_GoB] = orthogonalGoB_PrecodeGen(MV); % normalized precoder
R = nr_tot_SSBbeams*4;
tot_points = min(tot_points,R); % Q is set to be less than R
[beta_j_Pq, beta_j_Pq_fullloss] = beta_points(J,tot_points,point_positions,AP_positions,lambda); % fullloss, non-fullloss is the same
[beta_j_ue] = beta_j_ue_all(r,AP_positions,lambda,case_interf); % row=AP, col=links

% sensing channel
[h_ap1_2_pi, h_ap2_2_pi, h_ap3_2_pi, h_pi_2_aprx, h_APs_2_pi, h_direct_dep, h_direct_ari] = Channel_points2APs(MV,MH, ... 
    azi_deg_tx2points,ele_deg_tx2points,azi_deg_points2rx,ele_deg_points2rx, azi_direct_dep, azi_direct_ari);
% H1,.. where [M x M] % H_all=[H1;H2;H3] size is [M x JM x tot_points] % fullloss in evaluation
[H_all,H_all_fullloss,H01,H02,H03] = Gen_Hi(M,tot_points,h_pi_2_aprx,h_ap1_2_pi,h_ap2_2_pi,h_ap3_2_pi,beta_j_Pq, ... 
    beta_j_Pq_fullloss,radius_tx2points,radius_points2rx,fc,h_direct_dep, h_direct_ari,AP_positions,lambda); 
% combining vector v
V = 1/sqrt(M).*h_pi_2_aprx; %  v is [M x tot_points]

%list of give params % beta_j_ue
beta_clutter = 10^(-90/10); % -70 to -80 db
%% function for non-cordinated beamformer
h_direct_dep_n=h_direct_dep./norm(h_direct_dep); %normalized direct channel
rho_ssb_each = noise_var*gamma_ue/( beta_j_ue(1,1)- gamma_ue*sum(beta_j_ue(1,2:3)) ) ; % solve for rho_j
uncoordi_prec_AP1 = zeros(M,tot_points);
uncoordi_prec_AP2 = zeros(M,tot_points);
uncoordi_prec_AP3 = zeros(M,tot_points);
uncoordi_prec_all = zeros(M,J,tot_points);
compens=0;
for i=1:tot_points
    U_AP1 = [conj(precoderSSB_GoB(:,i+compens)),conj(h_direct_dep_n(:,1))]; % remove U = [conj(h_direct_dep_n(:,i))]
    U_AP2 = [conj(precoderSSB_GoB(:,i+compens)),conj(h_direct_dep_n(:,2))];
    U_AP3 = [conj(precoderSSB_GoB(:,i+compens)),conj(h_direct_dep_n(:,3))];
    uncoordi_prec_AP1(:,i) = (eye(M) - U*U')*h_ap1_2_pi(:,i);
    uncoordi_prec_AP2(:,i) = (eye(M) - U*U')*h_ap2_2_pi(:,i);
    uncoordi_prec_AP3(:,i) = (eye(M) - U*U')*h_ap3_2_pi(:,i);
    phase_change_AP2=exp(1i*2*pi*abs(radius_tx2points(2,i)-radius_tx2points(1,i)/lambda));
    phase_change_AP3=exp(1i*2*pi*abs(radius_tx2points(3,i)-radius_tx2points(1,i)/lambda));
    uncoordi_prec_AP1(:,i) = sqrt(Pmax-rho_ssb_each).*uncoordi_prec_AP1(:,i)./norm(uncoordi_prec_AP1(:,i));
    uncoordi_prec_AP2(:,i) = sqrt(Pmax-rho_ssb_each).*uncoordi_prec_AP2(:,i)./norm(uncoordi_prec_AP2(:,i)).*phase_change_AP2;
    uncoordi_prec_AP3(:,i) = sqrt(Pmax-rho_ssb_each).*uncoordi_prec_AP3(:,i)./norm(uncoordi_prec_AP3(:,i)).*phase_change_AP3;
    uncoordi_prec_all(:,1,i) = uncoordi_prec_AP1(:,i);
    uncoordi_prec_all(:,2,i) = uncoordi_prec_AP2(:,i);
    uncoordi_prec_all(:,3,i) = uncoordi_prec_AP3(:,i);
end

%% Optimization (CVX) for all points (len_point): Output=X_all [JM JM Q] rho_ssb_all [J Q] 
slack_t_db = 2; ebsilon = 0.2; % speed is one feasible is OK
[X_all,rho_ssb_all] = bisection_feasCheck_speed(tot_points,J,M,slack_t_db,ebsilon,H_all,V,precoderSSB_GoB,rcs_var, ... 
    beta_clutter,beta_j_ue,noise_var,gamma_ue,Pmax,h_direct_dep);

%% Conversion: Optimized X_all to precoder_vec [M J Q] ; s_pwr_mat [J Q] ; combine_s_pwr_p_ssb [J Q] (close to Pmax)
%gamma_ue=10^(3/10);
rho_ssb_each = noise_var*gamma_ue/( beta_j_ue(1,1)- gamma_ue*sum(beta_j_ue(1,2:3)) ) ; % solve for rho_j
rho_ssb_all=rho_ssb_each.*ones(J,tot_points);
[precoder_vec, s_pwr_mat, combine_s_pwr_p_ssb] = X_2_precoder_pwr_scale(M,J,tot_points,X_all,rho_ssb_all,Pmax);
% precoder_vec is scaled with power
%% EVALUATION of precoder in SINR wiht arbitary RCS
MC=1000; %len_point=1;
rcs_var_sim = 0.1; % arbitary RCS % MODIFY ME
rcs_dis = 1; % 1 SW2, 2 Weibull
[SICNR_sim_mean,SICNR_exp_mean,pwr_bgt_perAP_db,SICNR_sim_mean_uncoordi] = EvaluateSICNR(mod_order,MC,tot_points,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03,uncoordi_prec_all);

%% EVALUATION CDF
MC =1000; % MODIFY ME
rcs_dis=2; rcs_var_sim=0.1;
[sorted_SCINR_db, cdf_vals,precoder_vec_pwr] = EvaluateCDF_SINR(mod_order,MC,tot_points,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03);

%% EVALUATION Pd-Pfa, ROC
gamma_ue=10^(3/10);
rho_ssb_each = noise_var*gamma_ue/( beta_j_ue(1,1)- gamma_ue*sum(beta_j_ue(1,2:3)) ) ;% solve for rho_j
rho_ssb_all=rho_ssb_each.*ones(J,tot_points);
MC =3000; Peta = 0.1; %
rcs_dis=1; rcs_var_sim=0.1; pwr_bgt_perAP_db = 30; % Pmax 20,30 dB
[Pd_NC, Pfa_NC, Pd_WC, Pfa_WC,precoder_vec_pwr] = EvaluateROC(mod_order,MC,tot_points,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03,pwr_bgt_perAP_db,Peta);
Pd_WC
%% FIG4 ROC Pmax 21 dBm %%%%%%%%%%%%%%%%%%% 
figure
loglog(Pfa_NC,Pd_NC,'k-',LineWidth=2)
grid on
%legend(['$P' ...
%    '20 \; dBm$'], 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$P_{fa}$', 'Interpreter', 'latex', 'FontSize', 20),ylabel('$P_{d}$', 'Interpreter', 'latex', 'FontSize', 20)
hold on 

%% FIG4 ROC Pmax 20 dBm
loglog(Pfa_WC,Pd_WC,'r--',LineWidth=2)
%legend('$P_{max}=21 \; dBm, \; NC$','$P_{max}=20 \; dBm, \; WC$',{'$20 \; dBm$','$20 \; dBm$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%% FIG4 ROC Pmax 19 dBm
loglog(Pfa_WC,Pd_WC,'b-.',LineWidth=2)
%legend({'$P_{max}=20 \; dBm, \; NC$','$P_{max}=20 \; dBm, \; WC$','$P_{max}=30 \; dBm, \; NC$','$P_{max}=30 \; dBm, \; WC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%% FIG4 ROC Pmax 27 dBm NO COMBINER
loglog(Pfa_WC,Pd_WC,'m:',LineWidth=2)
%legend({'$P_{max}=21 \; dBm$','$P_{max}=20 \; dBm$','$P_{max}=19 \; dBm$','$P_{max}=27 \; dBm, \; NC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%% FIG4 ROC Pmax 21 dBm NO COMBINER
plot(Pfa_NC,Pd_NC,'k--',LineWidth=2) 
%legend({'$P_{max}=21 \; dBm$','$P_{max}=20 \; dBm$','$P_{max}=19 \; dBm$','$P_{max}=27 \; dBm, \; NC$','$P_{max}=21 \; dBm, \; NC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%% FIG4 ROC Pmax 21 dBm NO COMBINER
%plot(Pfa_NC,Pd_NC,'k--',LineWidth=2) 
%legend({'$P_{max}=26 \; dBm, \; WC, \; P_\eta=0.25$','$P_{max}=26 \; dBm, \; NC$','$P_{max}=20 \; dBm, \; WC$','$P_{max}=26 \; dBm, \; NC, \; P_\eta=0.25$','$P_{max}=20 \; dBm, \; NC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%legend({'$P_{max}=26 \; dBm, \; WC, \; P_\eta=0.25$','$P_{max}=26 \; dBm, \; NC$','$P_{max}=20 \; dBm, \; WC$','$P_{max}=26 \; dBm, \; NC, \; P_\eta=0.25$','$P_{max}=20 \; dBm, \; NC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 14)
legend({'$P_\mathrm{max}=30$ dBm $P_\eta=1$','$P_\mathrm{max}=40$ dBm $P_\eta=0.2$','$P_\mathrm{max}=30$ dBm $P_\eta=0.2$'}, 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 10)

%% FIG1 gamma 3,2db rcs 0.1 %%%%%%%%%%%%%%% z=10 in Fig1
figure()
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean ,'k-',LineWidth=2)
%hold on 
%semilogy(pwr_bgt_perAP_db,SCIR_exp_mean) % expectation
xlabel('$P_\mathrm{max}$ [dBm]','Interpreter', 'latex', 'FontSize', 16),ylabel('$\gamma_\mathrm{sen}$','Interpreter', 'latex', 'FontSize', 16)
grid on
%legend('\gamma = 3 dB, RCS: SW2, \sigma_{rcs}^2= 0.1', 'Location', 'northwest')
ylim([10^-3,10^5])
%% FIG1 gamma 15db rcs 0.1
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean_uncoordi,'b--',LineWidth=2)
%legend('\gamma = 3 dB, RCS:SW2, \sigma_{rcs}^2= 0.1','\gamma = 3 dB, RCS:SW2, \sigma_{rcs}^2= 0.001', 'Location', 'northwest')
%% FIG1 gamma 20 db rcs 0.1
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'r:',LineWidth=2)
%% Fig1
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean_uncoordi,'m-.',LineWidth=2)

xlabel('$P_\mathrm{max}$ (dBm)','Interpreter', 'latex', 'FontSize', 16),ylabel('$\gamma_\mathrm{sen}$','Interpreter', 'latex', 'FontSize', 18)
legend('$\gamma_{\mathrm{req}} = 3 $ dB, Proposed precoder', '$\gamma_{\mathrm{req}} = 3 $ dB, Non-coordinated precoder', ...
    '$\gamma_{\mathrm{req}} = 2 $ dB, Proposed precoder', '$\gamma_{\mathrm{req}} = 2 $ dB, Non-coordinated precoder', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12)

%% FIG2 gamma 3 db rcs 0.1 (z = 10, direct supressed) %%%%%%%%%%%%%%%%
figure()
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean ,'k-',LineWidth=2) %,Marker='o'
%hold on 
%semilogy(pwr_bgt_perAP_db,SCIR_exp_mean) % expectation
xlabel('$P_\mathrm{max}$ (dBm)','Interpreter', 'latex', 'FontSize', 16),ylabel('$\gamma_\mathrm{sen}$','Interpreter', 'latex', 'FontSize', 18)
%grid on
%% FIG2 gamma 3 db rcs 0.1 (z = 10, No direct supressed)
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'b--',LineWidth=2) % ,Marker='+'
%legend('\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.1','\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.001', 'Location', 'northwest')

%% FIG2 gamma 3 db rcs 0.1 (z = 1, direct supressed)
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'r:',LineWidth=2) % ,Marker='diamond'
%legend('$\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$','$\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.001$','$\gamma = 20 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12)

%% FIG2 gamma 3 db rcs 0.1 (z = 1, No direct supressed)
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'m-.',LineWidth=2) %,Marker='*'
legend('$z=1, $ DL-suppressed','$z=10, $ DL-suppressed', ... 
    '$z=10, $ NDL-suppressed' , ... 
    '$z=1, $ NDL-suppressed', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12)

xlim([18,25])
ylim([10^-1,10^3])

%% FIG 3 gamma 15db rcs 0.1 CDF
figure;
hold on
plot(sorted_SCINR_db, cdf_vals, 'k-', 'LineWidth', 2);
xlabel('$\gamma_{\mathrm{sen}}$','Interpreter', 'latex', 'FontSize', 18);
ylabel('CDF','Interpreter', 'latex', 'FontSize', 16);

%% FIG 3 gamma 15db rcs 0.001 CDF
hold on
plot(sorted_SCINR_db, cdf_vals, 'b--', 'LineWidth', 2);
legend('$\mathrm{SW2}, \sigma_{\mathrm{rcs}}^2= -10 \textrm{ dBsm}$' , ...
    '$\mathrm{WB}, w_\lambda=0.0847, w_k=0.833$', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12)

%% test plot gain pattern of the precoder
[vec, val] = eig(X_all(:,:,5)); 
mat_out_1 = conj(vec(:,end));
plotflag=1; edge_polar = pi/2;
[X,Y,Z,Gain_mat_abs] = plotBeampattern(mat_out_1(1:M),edge_polar,plotflag);
title('AP1')
%plotflag=1; edge_polar = pi/2;
[X,Y,Z,Gain_mat_abs] = plotBeampattern(mat_out_1(M+1:2*M),edge_polar,plotflag);
title('AP2')
%plotflag=1; edge_polar = pi/2;
[X,Y,Z,Gain_mat_abs] = plotBeampattern(mat_out_1(2*M+1:3*M),edge_polar,plotflag);
title('AP3')

%%

%function [Pd_sorted_WC,Pfa_sorted_WC,precoder_vec_pwr] = EvaluateROC(mod_order,MC,len_point,rho_ssb_all,rcs_var_sim, ... 
  %  precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03,pwr_bgt_perAP_db,Peta)
gamma_ue=10^(3/10);
rho_ssb_each = noise_var*gamma_ue/( beta_j_ue(1,1)- gamma_ue*sum(beta_j_ue(1,2:3)) ) ;% solve for rho_j
rho_ssb_all=rho_ssb_each.*ones(J,tot_points);
MC = 1e4; Peta = 0.05; %4.375e-4; %
rcs_var_sim=0.1; rcs_dis=1; pwr_bgt_perAP_db = 30; % Pmax 20,30 dB

SSB_symbols =  qammod(randi([0,mod_order-1],MC,J,len_point),mod_order,UnitAveragePower=true); % [MC Q] can extend to [MC J Q]
Sens_symbols =  qammod(randi([0,mod_order-1],MC,len_point),mod_order,UnitAveragePower=true); % [MC Q]

pwr_bgt_perAP = 10.^(pwr_bgt_perAP_db/10);
Sens_power_budget_mat = pwr_bgt_perAP - rho_ssb_all; %

T1_WC = zeros(len_point,MC);    % test statistics under H1
T0_WC = zeros(len_point,MC);    % test statistics under H0

%for mc=1:MC % MC rounds at each point
for q=1:len_point % Q iterate point
    q
    %T0_WC_sum=0;
    %T1_WC_sum=0;
    precoder_vec_pwr_AP1 = precoder_vec(:,1,q)./norm(precoder_vec(:,1,q)).*sqrt(Sens_power_budget_mat(1,q));
    precoder_vec_pwr_AP2 = precoder_vec(:,2,q)./norm(precoder_vec(:,2,q)).*sqrt(Sens_power_budget_mat(2,q));
    precoder_vec_pwr_AP3 = precoder_vec(:,3,q)./norm(precoder_vec(:,3,q)).*sqrt(Sens_power_budget_mat(3,q));
    %for q=1:len_point % Q iterate point
    for mc=1:MC % MC rounds at each point
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
                    ( H01*(Sens_symbols(mc,q).*precoder_vec_pwr_AP1 + SSB_symbols(mc,1,q).*sqrt(rho_ssb_all(1,q)).*conj(precoderSSB_GoB(:,q))) + ... %% direct link % both from SSB and Sensing (Sensing dominantate..)
                    H02*(Sens_symbols(mc,q).*precoder_vec_pwr_AP2 + SSB_symbols(mc,2,q).*sqrt(rho_ssb_all(2,q)).*conj(precoderSSB_GoB(:,q)) ) + ... 
                    H03*(Sens_symbols(mc,q).*precoder_vec_pwr_AP3 + SSB_symbols(mc,3,q).*sqrt(rho_ssb_all(3,q)).*conj(precoderSSB_GoB(:,q)) ) ); 
                eta_n = rand(1)<=Peta;
                y_received_T1 = eta_n.*Sens_symbols(mc,q).*( ... % sensing echoes [M x 1]
                    rcs1.*H_all_fullloss(:,1:M,q)*precoder_vec_pwr_AP1 + ... 
                    rcs2.*H_all_fullloss(:,M+1:2*M,q)*precoder_vec_pwr_AP2 + ... 
                    rcs3.*H_all_fullloss(:,2*M+1:3*M,q)*precoder_vec_pwr_AP3) + ...
                    y_received_T0;
                % combining V
                    y_combined_T0_WC = V(:,q)'*(y_received_T0);
                    y_combined_T1_WC = V(:,q)'*(y_received_T1);
                %end
                
                pwr_T0_WC = norm(y_combined_T0_WC)^2;
                pwr_T1_WC = norm(y_combined_T1_WC)^2; 

                %T0_WC_sum = T0_WC_sum+pwr_T0_WC;
                %T1_WC_sum = T1_WC_sum+pwr_T1_WC; 

                T0_WC(q,mc) = pwr_T0_WC;
                T1_WC(q,mc) = pwr_T1_WC; 
    end % points
    %T0_WC_sum_arr(mc)=T0_WC_sum;
    %T1_WC_sum_arr(mc)=T0_WC_sum;
end % mc

% ---- ROC computation ----
% thresholds spanning both distributions
Tmin_WC = min(min([T0_WC; T1_WC]));
Tmax_WC = max(max([T0_WC; T1_WC]));
%Tmin_WC_sum=min([T0_WC_sum_arr; T1_WC_sum_arr]);
%Tmax_WC_sum=max([T0_WC_sum_arr; T1_WC_sum_arr]);
nTh  = 2000;
Pd_WC = zeros(1,nTh);
Pfa_WC = zeros(1,nTh);
%Pd_sum = zeros(1,nTh);
%Pfa_sum = zeros(1,nTh);
taus_WC = linspace(Tmin_WC, Tmax_WC, nTh);
%taus_WC_sum = linspace(Tmin_WC_sum, Tmax_WC_sum, nTh);

for k = 1:nTh
    Pd_arr_WC = sum(T1_WC > taus_WC(k));
    Pd_WC(k)  = mean(Pd_arr_WC > 0);  % P_D = Pr(T > tau | H1) % one detect the whole sweep is detected
    Pfa_arr_WC = sum(T0_WC > taus_WC(k));
    Pfa_WC(k) = mean(Pfa_arr_WC > 0);  % P_FA = Pr(T > tau | H0)
    
    %Pd_sum(k) = mean(T0_WC_sum_arr> taus_WC_sum(k));
    %Pfa_sum(k) = mean(T1_WC_sum_arr> taus_WC_sum(k));
end

[Pfa_sorted_WC, idx_WC] = sort(Pfa_WC, 'ascend');
Pd_sorted_WC = Pd_WC(idx_WC);

%[Pfa_sorted_sum, idx_sum] = sort(Pfa_sum, 'ascend');
%Pd_sorted_sum = Pd_sum(idx_sum);

%% FIG4 ROC Pmax 21 dBm %%%%%%%%%%%%%%%%%%% 
figure
loglog(Pfa_sorted_WC,Pd_sorted_WC,'k-',LineWidth=2)
grid on
%legend(['$P' ...
%    '20 \; dBm$'], 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$P_{fa}$', 'Interpreter', 'latex', 'FontSize', 20),ylabel('$P_{d}$', 'Interpreter', 'latex', 'FontSize', 20)
%hold on 

%loglog(Pfa_sorted_sum,Pd_sorted_sum,'b-',LineWidth=2)

%% 
hold on
loglog(Pfa_sorted_WC,Pd_sorted_WC,'b--',LineWidth=2)
%% 
hold on
loglog(Pfa_sorted_WC,Pd_sorted_WC,'r:',LineWidth=2)
%% 
hold on
loglog(Pfa_sorted_WC,Pd_sorted_WC,'m-.',LineWidth=2)
%%
legend('$P_{\eta}=1$','$P_{\eta}=0.2$','$P_{\eta}=0.1$','$P_{\eta}=0.05$', 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 14)