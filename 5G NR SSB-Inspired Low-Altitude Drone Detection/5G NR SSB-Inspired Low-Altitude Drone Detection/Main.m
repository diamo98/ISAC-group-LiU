%% Initial Setup %%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%% setting: %%%%%%%%%%%%%
rcs_var_db = -10; gamma_ue_db = 15 ; z=10; Pmax = 10^(30/10); %%%%% MODIFY PARAMETERS here!
MV=12; MH=12; M=MV*MH; J=3; plot_flag = 1;
% %%%%% r-radius z-level, d-spacing, r_outer-hexagonal
r = 50; d=2; rcs_var = 10^(rcs_var_db/10); noise_var = 10^(-94/10);   
gamma_ue =  10^(gamma_ue_db/10); fc=15e9; lambda=3e8/fc; r_outer = r/cosd(30); % hexagonal: r is inner radius and r_outer is outer radius
% %%%%% define a space volumn [width x height x depth]
vol_width =2; vol_height = 2; vol_depth = 2;
point_row=ceil(vol_height/d +1); point_col=ceil(vol_width/d +1); point_deep=ceil(vol_depth/d +1);
num_points_tot = point_row*point_col*point_deep; Q = num_points_tot; mod_order = 4;

% define cartisian of points and APs 
[AP_positions,point_positions] = setup_cart(r,z,point_row,point_col,point_deep,d,plot_flag);
% define polar between APs_tx and points % azi_direct direct link
[azi_deg_tx2points, ele_deg_tx2points, radius_tx2points, azi_direct_dep, azi_direct_ari] = polar_APs_points(AP_positions,point_positions); 
% define polar between points and APs_rx
[azi_deg_points2rx, ele_deg_points2rx, radius_points2rx] = polar_Points_APrx(AP_positions,point_positions);

% W= [W1,W2,...,W_R] SSB beams and 
[azimu_GoB_SSB_deg, polar_GoB_SSB_deg, BeamList, nr_tot_SSBbeams, precoderSSB_GoB] = orthogonalGoB_PrecodeGen(MV); % normalized precoder
R = nr_tot_SSBbeams*4;
len_point = min(Q,R); % Q is set to be less than R
%precoder_Sens_APs = 1./sqrt(M).*conj(h_APs_2_pi); % [M x Q x 3] % to be optimized
%precoder_Sens_APs = ones(M,len_mat,3); % [M x Q x 3] % to be optimized

% beta_j_pq : the large scale fading
[beta_j_Pq, beta_SSB_echo, beta_j_Pq_fullloss] = beta_points(J,Q,point_positions,AP_positions,lambda);
%beta_j_ue =  [1./(4*pi*(x/2).^2).*ones(len_mat,1), 10^(-90/10).*ones(len_mat,2)];
[beta_j_ue] = beta_j_ue_all(r,AP_positions,lambda); % row=AP, col=links

%beta_j_Pq = ones(J,Q); % check precision issues
% define H=[H1,H2,..] at point-1 [MxQ]
[h_ap1_2_pi, h_ap2_2_pi, h_ap3_2_pi, h_pi_2_aprx, h_APs_2_pi, h_direct_dep, h_direct_ari] = Channel_points2APs(MV,MH, ... 
    azi_deg_tx2points,ele_deg_tx2points,azi_deg_points2rx,ele_deg_points2rx, azi_direct_dep, azi_direct_ari);
% H1, H2, H3 where [M x M x Q] % H_all=[H1;H2;H3] size is [M x JM x num_points_tot] % fullloss in evaluation
[H_all,H_all_fullloss,H01,H02,H03] = Gen_Hi(M,num_points_tot,h_pi_2_aprx,h_ap1_2_pi,h_ap2_2_pi,h_ap3_2_pi,beta_j_Pq, ... 
    beta_j_Pq_fullloss,radius_tx2points,radius_points2rx,fc,h_direct_dep, h_direct_ari,AP_positions,lambda); 
% combining vector v
V = 1/sqrt(M).*h_pi_2_aprx; %  v is [M x num_points_tot]

%list of give params % beta_j_ue
beta_clutter = 10^(-90/10); % -70 to -80 db
%mnoise = M*noise_var;
% first tier AP's rho to get SNR=10, gamma_ue_db ?? dB at cell edge of tier 1
beta_interfering_cell = lambda^2/(4*pi*(r_outer))^2;
tier1_rho = 10^(10/10)*noise_var/beta_interfering_cell;  % try first 
%tier1_rho=0;
%% Optimization (CVX) for all points (len_point): Output=X_all [JM JM Q] rho_ssb_all [J Q] 
slack_t_db = 2; ebsilon = 0.2; % speed is one feasible is OK
[X_all,rho_ssb_all] = bisection_feasCheck_speed(len_point,J,M,slack_t_db,ebsilon,H_all,V,precoderSSB_GoB,rcs_var, ... 
    beta_clutter,beta_j_ue,noise_var,gamma_ue,tier1_rho,Pmax,h_direct_dep);
% bisection_feasCheck

%% Conversion: Optimized X_all to precoder_vec [M J Q] ; s_pwr_mat [J Q] ; combine_s_pwr_p_ssb [J Q] (close to Pmax)
[precoder_vec, s_pwr_mat, combine_s_pwr_p_ssb] = X_2_precoder_pwr_scale(M,J,len_point,X_all,rho_ssb_all,Pmax);
% precoder_vec is scaled with power

%% EVALUATION of precoder in SINR wiht arbitary RCS
MC=3000; %len_point=1;
rcs_var_sim = 0.1; % arbitary RCS % MODIFY ME
rcs_dis = 1; % 1 SW2, 2 Weibull
%rho_ssb_all = 
[SICNR_sim_mean,SICNR_sim_mean_sep,SICNR_exp_mean,pwr_bgt_perAP_db] = EvaluateSICNR(mod_order,MC,len_point,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03);

%% EVALUATION CDF
MC =1000; % MODIFY ME
rcs_dis=2; rcs_var_sim=0.1; len_point=8;
[sorted_SCINR_db, cdf_vals] = EvaluateCDF_SINR(mod_order,MC,len_point,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03);

%% EVALUATION Pd-Pfa, ROC
MC =3000; Peta = 0.25; 
rcs_dis=1; rcs_var_sim=0.1; pwr_bgt_perAP_db = 26; % Pmax 20,30 dB
[Pd_NC, Pfa_NC, Pd_WC, Pfa_WC] = EvaluateROC(mod_order,MC,len_point,rho_ssb_all,rcs_var_sim, ... 
    precoder_vec,rcs_dis,J,M,H_all_fullloss,V,beta_clutter,precoderSSB_GoB,noise_var,H01,H02,H03,pwr_bgt_perAP_db,Peta);
Pd
%% FIG4 ROC Pmax 21 dBm %%%%%%%%%%%%%%%%%%% 
figure
plot(Pfa_WC,Pd_WC,'k--',Marker='+',LineWidth=2)
grid on
%legend(['$P' ...
%    '20 \; dBm$'], 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$P_{fa}$', 'Interpreter', 'latex', 'FontSize', 20),ylabel('$P_{d}$', 'Interpreter', 'latex', 'FontSize', 20)
hold on 
%% FIG4 ROC Pmax 20 dBm
plot(Pfa_NC,Pd_NC,'m--',Marker='*',LineWidth=2)
%legend('$P_{max}=21 \; dBm, \; NC$','$P_{max}=20 \; dBm, \; WC$',{'$20 \; dBm$','$20 \; dBm$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%% FIG4 ROC Pmax 19 dBm
plot(Pfa_WC,Pd_WC,'b-.',Marker='o',LineWidth=2)
%legend({'$P_{max}=20 \; dBm, \; NC$','$P_{max}=20 \; dBm, \; WC$','$P_{max}=30 \; dBm, \; NC$','$P_{max}=30 \; dBm, \; WC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%% FIG4 ROC Pmax 27 dBm NO COMBINER
plot(Pfa_NC,Pd_NC,'r:',Marker='x',LineWidth=2)
%legend({'$P_{max}=21 \; dBm$','$P_{max}=20 \; dBm$','$P_{max}=19 \; dBm$','$P_{max}=27 \; dBm, \; NC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%% FIG4 ROC Pmax 21 dBm NO COMBINER
plot(Pfa_NC,Pd_NC,'k--',LineWidth=2) 
%legend({'$P_{max}=21 \; dBm$','$P_{max}=20 \; dBm$','$P_{max}=19 \; dBm$','$P_{max}=27 \; dBm, \; NC$','$P_{max}=21 \; dBm, \; NC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
%% FIG4 ROC Pmax 21 dBm NO COMBINER
%plot(Pfa_NC,Pd_NC,'k--',LineWidth=2) 
%legend({'$P_{max}=26 \; dBm, \; WC, \; P_\eta=0.25$','$P_{max}=26 \; dBm, \; NC$','$P_{max}=20 \; dBm, \; WC$','$P_{max}=26 \; dBm, \; NC, \; P_\eta=0.25$','$P_{max}=20 \; dBm, \; NC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
legend({'$P_{max}=26 \; dBm, \; WC, \; P_\eta=0.25$','$P_{max}=26 \; dBm, \; NC$','$P_{max}=20 \; dBm, \; WC$','$P_{max}=26 \; dBm, \; NC, \; P_\eta=0.25$','$P_{max}=20 \; dBm, \; NC$'}, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 14)

%% FIG1 gamma 15db rcs 0.1 %%%%%%%%%%%%%%% z=10 in Fig1
figure()
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean ,'k-',LineWidth=2)
%hold on 
%semilogy(pwr_bgt_perAP_db,SCIR_exp_mean) % expectation
xlabel('Pmax [dBm]','Interpreter', 'latex', 'FontSize', 16),ylabel('SINR','Interpreter', 'latex', 'FontSize', 16)
grid on
legend('\gamma = 15 dB, RCS: SW2, \sigma_{rcs}^2= 0.1', 'Location', 'northwest')
%% FIG1 gamma 15db rcs 0.001
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'b--',LineWidth=2)
legend('\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.1','\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.001', 'Location', 'northwest')
%% FIG1 gamma 20 db rcs 0.1
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'r:',LineWidth=2)
legend('$\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$','$\gamma = 20 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$', ...
    '$\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.001$' , 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('Pmax [dBm]','Interpreter', 'latex', 'FontSize', 18),ylabel('SINR','Interpreter', 'latex', 'FontSize', 18)
%% FIG1 gamma 20 db weibull weibull_lambda = 0.0847;weibull_k = 0.833;
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'m-.',LineWidth=2)
legend('$\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$','$\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.001$', ...
    '$\gamma = 20 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$' , ...
    '$\gamma = 20 dB, RCS:WB, w_\lambda=0.0847, w_k=0.833$', 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Pmax [dBm]','Interpreter', 'latex', 'FontSize', 20),ylabel('SINR','Interpreter', 'latex', 'FontSize', 20)



%% FIG2 gamma 15 db rcs 0.1 (z = 10, direct supressed) %%%%%%%%%%%%%%%%
figure()
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean ,'k-',LineWidth=2) %,Marker='o'
%hold on 
%semilogy(pwr_bgt_perAP_db,SCIR_exp_mean) % expectation
xlabel('Pmax [dBm]','Interpreter', 'latex', 'FontSize', 20),ylabel('SINR','Interpreter', 'latex', 'FontSize', 20)
grid on
%% FIG2 gamma 15 db rcs 0.1 (z = 10, No direct supressed)
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'r--',LineWidth=2) % ,Marker='+'
%legend('\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.1','\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.001', 'Location', 'northwest')

%% FIG2 gamma 15 db rcs 0.1 (z = 1, direct supressed)
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'b:',LineWidth=2) % ,Marker='diamond'
%legend('$\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$','$\gamma = 15 dB, RCS:SW2, \sigma_{rcs}^2= 0.001$','$\gamma = 20 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12)

%% FIG2 gamma 15 db rcs 0.1 (z = 1, No direct supressed)
hold on
semilogy(pwr_bgt_perAP_db,SICNR_sim_mean,'k:.',Marker='+',LineWidth=2) %,Marker='*'
legend('$z=10, DL$-$suppressed$','$z=1, DL$-$suppressed$', ... 
    '$z=10, N$-$DL$-$suppressed$' , ... 
    '$z=1, N$-$DL$-$suppressed$', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12)

%% FIG 3 gamma 15db rcs 0.1 CDF
figure;
hold on
plot(sorted_SCINR_db, cdf_vals, 'k-', 'LineWidth', 2);
xlabel('SINR [dB]','Interpreter', 'latex', 'FontSize', 16);
ylabel('CDF','Interpreter', 'latex', 'FontSize', 16);
grid on;
%% FIG 3 gamma 15db rcs 0.001 CDF
hold on
plot(sorted_SCINR_db, cdf_vals, 'b--', 'LineWidth', 2);
%% FIG 3 gamma 20 db rcs 0.001 CDF

plot(sorted_SCINR_db, cdf_vals, 'r:', 'LineWidth', 2);
%% FIG 3 gamma 20 db rcs WB CDF
hold on
plot(sorted_SCINR_db, cdf_vals, 'm-.', 'LineWidth', 2);
grid on
xlabel('SINR [dB]','Interpreter', 'latex', 'FontSize', 20), ylabel('CDF','Interpreter', 'latex', 'FontSize', 20);
legend('$\gamma = 20 dB, RCS:SW2, \sigma_{rcs}^2= 0.1$' , ...
    '$\gamma = 20 dB, RCS:WB, w_\lambda=0.0847, w_k=0.833$', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 16)

%% test plot gain pattern of the precoder
[vec, val] = eig(X_all(:,:,1)); 
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

