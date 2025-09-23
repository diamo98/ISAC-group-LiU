%% SETTING 
M=100; % ISAC MIMO antennas
N=100; % number of repeaters
theta_drone = pi/6; % drone angle rad
theta_drone_deg = rad2deg(theta_drone); % drone angle deg
lambda = 3e8/15e9; rcs= 10^(-10/10); 
% function
a = @(theta) [exp(1i*pi*(0:M-1)'*cos(theta))]; % array response
onelink_beta = @(l1) (lambda^2./(4*pi.*l1).^2);
twolink_beta = @(l1,l2) (rcs*lambda^2./((4*pi)^3.*l1.^2.*l2.^2));
% position
frac_dis = 1; % 1-regular range, 10-short range
l_AD = 500/frac_dis ; l_A1 = 250/frac_dis ; l_AU = 100/frac_dis ; d = 400/N/frac_dis ; % repeater spacing
AP_pos = [0;0]; D_pos = [l_AD*cos(theta_drone) ; l_AD*sin(theta_drone)]; Repeater_pos = [l_A1+(0:N-1)*d;zeros(1,N)];
l_AnR = l_A1+(0:N-1)*d; l_DnR = vecnorm(Repeater_pos-D_pos);
% beta
beta_AD= twolink_beta(l_AD,l_AD); beta_AU= onelink_beta(l_AU); beta_AnR= onelink_beta(l_AnR)'; beta_ADnR = twolink_beta(l_AD,l_DnR)';
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

% %%%%%%% TX signal
MC=1000; mod=4;
com_sym=randQPSK(MC);
sen_sym=randQPSK(MC);
W_s_p=conj(a(theta_drone))/norm(a(theta_drone));
[rho_c, rho_s] = optimal_power_Split(gamma_UE,M,sigma_ue2,beta_AU,rho_max);
% channel
delay = @(dis) exp(-1i*2*pi*dis/3e8);
H_ad = sqrt(beta_AD)*delay(l_AD*2)*a(theta_drone)*conj(a(theta_drone)');
H_ar = (a(0).*ones(1,N))*diag(delay(l_AnR).*sqrt(beta_AnR)');
H_adr = (a(theta_drone).*ones(1,N)) *diag(delay(l_AD+l_DnR).*sqrt(beta_ADnR)');
H_rr = 0;
phi_alpha = diag(alpha);


% for i=1:MC % tx signal
%     h_au = sqrt(beta_AU/2)*(randn(M,1)+1i*randn(M,1)); % random channel of a user 
%     W_s= (eye(M)-h_au*h_au')*W_s_p; 
%     n_r = sqrt(sigma_r2/2)*randn(N,1)+1i*randn(N,1); % repeater noise
%     n_ap_ue = sqrt(sigma_ap2_ue2/2)*randn(M,1)+1i*randn(M,1); % AP noise
%     W_c= conj(h_au)/norm(h_au); % communication precoder
%     X = sqrt(rho_c)*W_c*com_sym(i)+sqrt(rho_s)*W_s*sen_sym(i);
%     Y = H_ad*X + H_ar*phi_alpha*conj(H_adr')*X + H_ar*phi_alpha*(conj(H_adr')*X + n_r) + n_ap_ue;
% 
% end

% ----- Inputs you already have (example placeholders) -----
% M, N, MC, theta, beta_AU, sigma_r2, sigma_ap2_ue2, rho_c, rho_s,
% H_ad, H_ar, H_adr, phi_alpha (N x N diag of alphas), W_s_p
% com_sym (MC x 1), sen_sym (MC x 1)  -- if not, we'll generate QPSK.

% Helper: ULA steering a(theta) = [1, e^{jπcosθ}, ..., e^{j(M-1)πcosθ}]^T
steer = @(th,M) exp(1j*pi*cos(th)*(0:M-1)).';

%% ROC for each line with setting [alpha_max,N]
MC = 5000;
if ~exist('com_sym','var') || numel(com_sym) ~= MC
    b = randi([0 1], MC, 2);  % QPSK
    com_sym = ((1-2*b(:,1)) + 1j*(1-2*b(:,2))) / sqrt(2);
end
if ~exist('sen_sym','var') || numel(sen_sym) ~= MC
    b = randi([0 1], MC, 2);
    sen_sym = ((1-2*b(:,1)) + 1j*(1-2*b(:,2))) / sqrt(2);
end

% Matched filter vector toward target AoA
a_theta = steer(theta_drone, M);
w = conj(a_theta) / norm(a_theta);   % unit-norm beam

% Preallocate
Y1 = zeros(M, MC);   % H1: target present
Y0 = zeros(M, MC);   % H0: target absent
T1 = zeros(MC,1);    % test statistics under H1
T0 = zeros(MC,1);    % test statistics under H0

% Zeroed target channels for H0
H_ad0  = zeros(M,M);
H_adr0 = zeros(M,N);

% ---- Monte Carlo ----
for i = 1:MC
    % user channel
    h_au = sqrt(beta_AU/2) * (randn(M,1) + 1j*randn(M,1));

    % sensing precoder orthogonal to h_au
    W_s = (eye(M) - h_au*h_au'/(norm(h_au)^2 + eps)) * W_s_p;
    norm(W_s)
    % noises
    n_r    = sqrt(sigma_r2/2)      * (randn(N,1) + 1j*randn(N,1));   % repeater noise
    n_ap   = sqrt(sigma_ap2_ue2/2) * (randn(M,1) + 1j*randn(M,1));   % AP noise

    % comm precoder
    W_c = conj(h_au) / (norm(h_au) + eps);

    % transmit vector
    X = sqrt(rho_c)*W_c*com_sym(i) + sqrt(rho_s)*W_s*sen_sym(i);

    % ----- H1: target present (your original model) -----
    Y = H_ad*X ...
      + H_ar*phi_alpha*conj(H_adr')*X ...
      + H_ar*phi_alpha*(conj(H_adr')*X + n_r) ...
      + n_ap;
    Y1(:,i) = Y;

    % ----- H0: target absent (remove drone reflections) -----
    Y_0 = H_ad0*X ...
        + H_ar*phi_alpha*conj(H_adr0')*X ...
        + H_ar*phi_alpha*(conj(H_adr0')*X + n_r) ...
        + n_ap;
    Y0(:,i) = Y_0;

    % Matched-filter test statistic
    T1(i) = abs(w' * Y)^2;
    T0(i) = abs(w' * Y_0)^2;
end

% ---- ROC computation ----
% thresholds spanning both distributions
Tmin = min([T0; T1]);
Tmax = max([T0; T1]);
nTh  = 400;
taus = linspace(Tmin, Tmax, nTh);

Pd  = zeros(nTh,1);
Pfa = zeros(nTh,1);
for k = 1:nTh
    tau = taus(k);
    Pd(k)  = mean(T1 > tau);  % P_D = Pr(T > tau | H1)
    Pfa(k) = mean(T0 > tau);  % P_FA = Pr(T > tau | H0)
end

% Optional: area under curve
% Sort Pfa ascending and reorder Pd accordingly
[Pfa_sorted, idx] = sort(Pfa(:), 'ascend');
%Pd_sorted = Pd(idx);
AUC = trapz(sort(Pfa), Pd(idx));

%% ---- Plot ROC for each case [alpha_max,N]----
%plot(Pfa, Pd,'k-', 'LineWidth', 2); hold on; % N=100
plot(Pfa, Pd,'k--', 'LineWidth', 2); hold on; % N=10
%plot(Pfa, Pd,'r:', 'LineWidth', 2); hold on; % N=0 alpha=40
%plot(Pfa, Pd,'r-.', 'LineWidth', 2); hold on; % N=10
%plot(Pfa, Pd,'b:', 'LineWidth', 2); hold on; % N=0 alpha=40

grid on; xlim([0 1]); ylim([0 1]);

hold on

%% Labels and legend
xlabel('$P_{FA}$','FontSize',18,'Interpreter','latex');
ylabel('$P_{D}$','FontSize',18,'Interpreter','latex');
%plot([0 1],[0 1],'k--','LineWidth',1);
legend('N=100, $\alpha_{max}=48$ dB','N=100, $\alpha_{max}=45$ dB','N=50, $\alpha_{max}=48$ dB','N=50, $\alpha_{max}=45$ dB','$N=0$ dB','FontSize',14,'Interpreter','latex','Location', 'southeast')
% ---- (Optional) CFAR threshold for a desired P_FA0 ----
% P_FA0 = 0.05;  % target false alarm
% tau_cfar = quantile(T0, 1 - P_FA0);
% Pd_cfar  = mean(T1 > tau_cfar);
% fprintf('CFAR threshold for P_FA=%.3f -> tau=%.3g, Pd=%.3f\n', P_FA0, tau_cfar, Pd_cfar);
% 






