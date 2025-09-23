% ISAC NCR swarm
clear all
M=100; % ISAC MIMO antennas
N=50; % number of repeaters
theta_drone = pi/6; % drone angle rad
theta_drone_deg = rad2deg(theta_drone); % drone angle deg
lambda = 3e8/15e9; rcs= 10^(-10/10); 
% function
a = @(theta) [exp(1i*pi*(0:M-1)'*cos(theta))]; % array response
onelink_beta = @(l1) (lambda^2./(4*pi.*l1).^2);
twolink_beta = @(l1,l2) (rcs*lambda^2./((4*pi)^3.*l1.^2.*l2.^2)) ;
% position
frac_dis = 1;
l_AD = 500/frac_dis ; l_A1 = 250/frac_dis ; l_AU = 100/frac_dis ; d = 400/N/frac_dis ; % repeater spacing
AP_pos = [0;0]; D_pos = [l_AD*cos(theta_drone) ; l_AD*sin(theta_drone)]; Repeater_pos = [l_A1+(0:N-1)*d;zeros(1,N)];
l_AnR = l_A1+(0:N-1)*d; l_DnR = vecnorm(Repeater_pos-D_pos);
% beta
beta_AD = twolink_beta(l_AD,l_AD); beta_AU= onelink_beta(l_AU); beta_AnR= onelink_beta(l_AnR)'; beta_ADnR = twolink_beta(l_AD,l_DnR)';
% noise
sigma_r2 = 10^(-124/10); sigma_ap2_ue2 = 10^(-110/10); sigma_ue2 = sigma_ap2_ue2; sigma_ap2 = sigma_ap2_ue2;
gamma_UE = 10^(15/10); % UE_SINR threshold 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_max = 10^(45/10); rho_max = 10^(33/10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[alpha, rho_c, rho_s, gamma_opt, info] = solve_NCR_gain( ...
  beta_AD, beta_AU, beta_AnR, beta_ADnR, ...
  sigma_r2, sigma_ap2_ue2, sigma_ap2_ue2, gamma_UE, ...
  M, rho_max, alpha_max, struct('verbose',1));

alpha'

%% Evaluate Sensing SINR over Pmax (rho_max) [not used]
% Weak Sensing Channel, alpha_max [0,50,75,90]
rho_max_arr_dBm = [-20:40];
rho_max_arr = 10.^(rho_max_arr_dBm./10);
alpha_max_arr_dB = [-inf,50,75,90];
alpha_max_arr = 10.^(alpha_max_arr_dB/10);
SINR_sens_arr = zeros(length(alpha_max_arr_dB),length(rho_max_arr));
for j=1:length(alpha_max_arr)
    
    for i=1:length(rho_max_arr)
        [alpha, rho_c, rho_s, gamma_opt, info] = solve_NCR_gain( beta_AD, beta_AU, beta_AnR, beta_ADnR, sigma_r2, sigma_ap2_ue2, sigma_ap2_ue2, gamma_UE, M, rho_max_arr(i), alpha_max_arr(j), struct('verbose',1));
        %[rho_c, rho_s] = optimal_power_Split(gamma_UE,M,sigma_ue2,beta_AU,rho_max_arr(i));
        gamma = sensing_gamma(alpha, rho_c, rho_s, M, ...
                                       beta_AD, beta_AnR, beta_ADnR, ...
                                       sigma_r2, sigma_ap2_ue2);
        SINR_sens_arr(j,i) = gamma;
    end
end


figure()
semilogy(rho_max_arr_dBm,SINR_sens_arr(1,:),'k-','LineWidth',1.5)
hold on
semilogy(rho_max_arr_dBm,SINR_sens_arr(2,:),'r--','LineWidth',1.5)
semilogy(rho_max_arr_dBm,SINR_sens_arr(3,:),'b:','LineWidth',1.5)
semilogy(rho_max_arr_dBm,SINR_sens_arr(4,:),'k--','LineWidth',1.5,"Marker","diamond") %,"Marker","diamond"
grid on
legend('$N=0$', '$N=50,\alpha_{\max}=50$', '$N=50,\alpha_{\max}=75$', '$N=50,\alpha_{\max}=90$', ...
       'Interpreter', 'latex', ...
       'FontSize', 12, ...
       'Location', 'southeast');

%% Evaluate Sensing SINR over Alpha_max WC
% Weak Sensing Channel, 
N_arr = [0,5,10,100];
alpha_max_arr_dB = [20:0.5:90];
alpha_max_arr = 10.^(alpha_max_arr_dB/10);

SINR_sens_arr = zeros(length(N_arr),length(alpha_max_arr_dB));
track_max_alpha = zeros(length(N_arr),max(N_arr));
frac_dis = 1;
for i=1:length(N_arr)
    
    for j=1:length(alpha_max_arr)
        % change setting
        [beta_AD, beta_AU, beta_AnR, beta_ADnR, ] = modify_scheme(N_arr(i),frac_dis) ;
        %[alpha, rho_c, rho_s, gamma_opt, info] = solve_NCR_gain(beta_AD, beta_AU, beta_AnR, beta_ADnR,sigma_r2, sigma_ap2_ue2, sigma_ap2_ue2, gamma_UE, M, rho_max, alpha_max, struct('verbose',1));
        [alpha, rho_c, rho_s, gamma_opt, info] = solve_NCR_gain( beta_AD, beta_AU, beta_AnR, beta_ADnR, sigma_r2, sigma_ap2_ue2, sigma_ap2_ue2, gamma_UE, M, rho_max, alpha_max_arr(j), struct('verbose',1));
        % [rho_c, rho_s] = optimal_power_Split(gamma_UE,M,sigma_ue2,beta_AU,rho_max_arr(i));
        % non_zero_arr(j)=nnz(alpha);
        track_max_alpha(i,1:N_arr(i)) = ((alpha'>track_max_alpha(i, 1:N_arr(i) )).*alpha')+ (track_max_alpha(i, 1:N_arr(i) ).*~(alpha'>track_max_alpha(i, 1:N_arr(i) ))); % this keeps the alpha max for each repeater n
        %track_max_alpha(i,1:N_arr(i))
        gamma = sensing_gamma(track_max_alpha(i,1:N_arr(i)) , rho_c, rho_s, M, ... 
                                       beta_AD, beta_AnR, beta_ADnR, ... 
                                       sigma_r2, sigma_ap2_ue2);

        SINR_sens_arr(i,j) = gamma;

    end
end

%
figure()
semilogy(alpha_max_arr_dB,SINR_sens_arr(4,:),'k-','LineWidth',1.8)
hold on
semilogy(alpha_max_arr_dB,SINR_sens_arr(3,:),'b--','LineWidth',1.8)
semilogy(alpha_max_arr_dB,SINR_sens_arr(2,:),'R:','LineWidth',1.8)
semilogy(alpha_max_arr_dB,SINR_sens_arr(1,:),'g:','LineWidth',1.8)
grid on
xlabel('$\alpha_{max} \; [dB]$','FontSize',18,'Interpreter','latex'),ylabel('$\gamma_s$','FontSize',18,'Interpreter','latex')
legend('$N=100$', '$N=10$', '$N=5$', '$N=0$' , 'Interpreter', 'latex', 'FontSize', 14,'Location', 'northwest');

%% Evaluate Sensing SINR over Alpha_max SC
% strong Sensing Channel, 
N_arr = [0,5,10,100];
alpha_max_arr_dB = [20:0.5:90];
alpha_max_arr = 10.^(alpha_max_arr_dB/10);

SINR_sens_arr = zeros(length(N_arr),length(alpha_max_arr_dB));
track_max_alpha = zeros(length(N_arr),max(N_arr));
frac_dis = 10;
for i=1:length(N_arr)
    
    for j=1:length(alpha_max_arr)
        % change setting
        [beta_AD, beta_AU, beta_AnR, beta_ADnR, ] = modify_scheme(N_arr(i),frac_dis) ;
        [alpha, rho_c, rho_s, gamma_opt, info] = solve_NCR_gain( beta_AD, beta_AU, beta_AnR, beta_ADnR, sigma_r2, sigma_ap2_ue2, sigma_ap2_ue2, gamma_UE, M, rho_max, alpha_max_arr(j), struct('verbose',1));
        % [rho_c, rho_s] = optimal_power_Split(gamma_UE,M,sigma_ue2,beta_AU,rho_max_arr(i));
        % non_zero_arr(j)=nnz(alpha);
        track_max_alpha(i,1:N_arr(i)) = ((alpha'>track_max_alpha(i, 1:N_arr(i) )).*alpha')+ (track_max_alpha(i, 1:N_arr(i) ).*~(alpha'>track_max_alpha(i, 1:N_arr(i) ))); % this keeps the alpha max for each repeater n

        gamma = sensing_gamma(track_max_alpha(i,1:N_arr(i)) , rho_c, rho_s, M, ... 
                                       beta_AD, beta_AnR, beta_ADnR, ... 
                                       sigma_r2, sigma_ap2_ue2)
        SINR_sens_arr(i,j) = gamma;

    end
end
%
figure()
semilogy(alpha_max_arr_dB,SINR_sens_arr(4,:),'k-','LineWidth',1.8)
hold on
semilogy(alpha_max_arr_dB,SINR_sens_arr(3,:),'b--','LineWidth',1.8)
semilogy(alpha_max_arr_dB,SINR_sens_arr(2,:),'r:','LineWidth',1.8)
semilogy(alpha_max_arr_dB,SINR_sens_arr(1,:),'g:','LineWidth',1.8)
grid on
xlabel('$\alpha_{max} \; [dB]$','FontSize',18,'Interpreter','latex'),ylabel('SINR ($\gamma_s$)','FontSize',18,'Interpreter','latex')
legend('$SC, N=100$', '$SC, N=10$', '$SC, N=5$', '$SC, N=0$','Interpreter', 'latex', 'FontSize', 14,'Location', 'southeast');

