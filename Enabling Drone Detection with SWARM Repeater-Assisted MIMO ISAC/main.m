%% Enabling Drone Detection with SWARM Repeater-Assisted MIMO ISAC
% SETUP %
clear all
addpath('./utils');
M=100; % ISAC MIMO antennas
%........................
N=10; % number of repeaters
d = 100 ; % repeater spacing 
%........................
theta_drone = pi/6; % drone angle rad 
c = 3e8;
theta_drone_deg = rad2deg(theta_drone); % drone angle in degree
fc = 15e9; lambda = c/fc; sigma_rcs2= 10^(-10/10);  % carrier frequency, wavelength, rcs variance
% array response of the target for perfect DFT beam (angle is matched)
a = @(theta) [exp(1i*pi*(0:M-1)'*cos(theta))]; 
onelink_beta = @(l1) (lambda^2./(4*pi.*l1).^2); % define function pathloss/channel gain
twolink_beta = @(l1,l2) (lambda^2./((4*pi)^3.*l1.^2.*l2.^2)) ; % bistatic link pathloss
% Reater positioning
frac_dis = 1; % fraction coeffcicient for repeater spacing, default is 1
l_AD = 500 ; % AP - Drone distance
l_A1 = 250/frac_dis ; % AP - first repeater distance
l_AU = 200/frac_dis ;  % AP - user distance
AP_pos = [0;0]; % AP position in 2D coordinate
D_pos = [l_AD*cos(theta_drone) ; l_AD*sin(theta_drone)]; % Drone position in 2D coordinate
Repeater_pos = [l_A1+(0:N-1)*d ; zeros(1,N)];  % Repeaters position
l_AnR = l_A1+(0:N-1)*d;  % array contains distance between AP - repeater n-th
l_DnR = vecnorm(Repeater_pos-D_pos); % array contains distace between Drone - repeater n-th
tau_AD = l_AD*2/c; % time propagation AP - Drone
tau_AnR = l_AnR/c; % time propagation AP - repeater n-th
tau_ADn = (l_AD+l_DnR)/c; % time propagation AP - Drone - repeater n-th
beta_AD = twolink_beta(l_AD,l_AD); % monostatic link pathloss
beta_AU = onelink_beta(l_AU);  % channel gain AP - User
beta_AR = onelink_beta(l_AnR)'; % channel gain AP - repeater n-th
beta_ADR = twolink_beta(l_AD,l_DnR)'; % bistatic link pathloss
sigma_R2 = 10^(-110/10); % repeater noise variance
sigma_AP2 = 10^(-90/10); % AP and UE noise variance
%% Power budget AP, User SIRN constraint
rho_max = 10^(53/10); % 200 watts, AP power budget
gamma_UE_req = 10^(30/10); % user SINR constraint
rho_c = gamma_UE_req*sigma_AP2/(M*beta_AU); % communication power
rho_s = rho_max - rho_c; % sensing power
%% define channels and alpha_max
w_s_p = a(theta_drone)./norm(a(theta_drone)); % normalized steering vector for target sensing
H_AD = sqrt(beta_AD).*exp(-1i*2*pi*fc*tau_AD).*a(theta_drone)*conj(a(theta_drone)'); % channel  AP-drone [M M]
H_AR = sqrt(beta_AR)'.*exp(-1i*2*pi*fc*tau_AnR).*a(0); % channel  AP-repeater [M N]
H_ADR = sqrt(beta_ADR)'.*exp(-1i*2*pi*fc*tau_ADn).*a(theta_drone); % chennel  AP-drone-repeater [M N]
H_RR = gen_H_RR(N,Repeater_pos); % inter-repeater channel, without self interference

alpha_max = min(1./sum(abs(H_RR))); % maximum amplification gain for stability system
alpha_max_dB = 10*log10(alpha_max); % in dB
fprintf('Alpha max is: %d , in dB : %d \n',alpha_max,alpha_max_dB)

%% This section validate empirical SINR and closed-form SINR
ALPHA = alpha_max*1/2*eye(N); % test value 
%ALPHA = 0*eye(N); % test value 
mc=5000; gamma_sen_sim_nu = 0;  gamma_sen_sim_denom = 0; temp1=0; temp2=0; temp3 = 0; temp4=0; temp5 =0;

for i=1:mc
    h_au = sqrt(beta_AU/2).*(randn(M,1)+1i*randn(M,1));  % randomize user channel form channel estimate
    A1 = conj(h_au)./norm(h_au); A2 = conj(a(0))./norm(a(0));
    U = [A1,A2];
    w_s = (eye(M) - U/(U'*U)*U')*conj(w_s_p) ; % sensing precoder 
    % w_s = conj(w_s_p); % do not use, perform veary bad in sensing SINR 
    w_s = (w_s)./norm(w_s); % normalized sensing precoder
    A0_n = conj(a(0))*conj(a(0))'./norm(a(0))^2 ;
    w_c = (eye(M) - A0_n)*conj(h_au); % communication precoder
    % w_c = conj(h_au); % do not use, perform veary bad in sensing SINR 
    w_c = (w_c)./norm(w_c); % normalzied communication precoder
    s_s = qammod(randi([0,3],1),4,'UnitAveragePower',true); s_c = qammod(randi([0,3],1),4,'UnitAveragePower',true); %s_s=1; s_c=1;
    x = sqrt(rho_s)*w_s*s_s + sqrt(rho_c)*w_c*s_c ;
    noise_R = sqrt(sigma_R2/2)*(randn(N,1)+1i*randn(N,1));
    noise_AP = sqrt(sigma_AP2/2)*(randn(M,1)+1i*randn(M,1));
    noise_UE = sqrt(sigma_AP2/2)*(randn()+1i*randn());
    rcs_reli = sqrt(sigma_rcs2/2)*(randn(1,1)+1i*randn(1,1));
    gamma_sen_sim_nu = gamma_sen_sim_nu + norm(rcs_reli.*H_AD*x + H_AR*(eye(N) + ALPHA*H_RR)*ALPHA*conj(H_ADR')*x.*rcs_reli)^2; % numerator 
    % gamma_sen_sim_nu = gamma_sen_sim_nu + norm(H_AR*(eye(N) + ALPHA*H_RR)*ALPHA*conj(H_ADR')*x)^2; % numerator 
    
    % gamma_sen_sim_denom = gamma_sen_sim_denom + norm(H_AR*ALPHA*(conj(H_AR)'*x + noise_R) + noise_AP)^2 ; 
    gamma_sen_sim_denom = gamma_sen_sim_denom + norm(H_AR*(eye(N) + ALPHA*H_RR)*ALPHA*noise_R + noise_AP)^2 ; 
    
    % sanity CHECK of SINR at each component
    %temp1 = temp1 + norm(H_AD*x)^2; 
    % noise_only = noise_only + norm(noise_AP)^2; 
    %temp2 = temp2 + norm((H_AR*ALPHA)*conj(H_ADR)'*x)^2; 
    %temp3 = temp3 + norm((H_AR/(eye(N)-ALPHA*H_RR)*ALPHA)*conj(H_ADR)'*x)^2;
    %temp4 = temp4 + abs(conj(h_au')*sqrt(rho_s)*w_s)^2;
    %temp5 = temp5 + abs(rcs_reli)^2;
end

empirical = gamma_sen_sim_nu/gamma_sen_sim_denom;
fprintf('empiricalsensing SINR is: %d \n',empirical)
% closed form Eq.19a
c1 = sigma_rcs2 * beta_AD * (rho_c + rho_s * M);
c2 = sigma_rcs2 * (rho_c + rho_s * M);
d_vec = sqrt(beta_AR);
v = sqrt(beta_ADR).*exp(-1i*2*pi*fc*tau_ADn'); % Repeater to Target channel
B = (eye(N) + ALPHA*H_RR)*ALPHA;
Numerator_closeform = c1 + c2*abs(conj(d_vec)'*B*v )^2;
Denomenator_closedform = sigma_R2*norm(B'*d_vec)^2 + sigma_AP2;
closeform_SINR = Numerator_closeform/Denomenator_closedform;
fprintf('Closed-form SINR is: %d \n',closeform_SINR)

%% Optimization ALPHAs WITH different rho_max
rho_max_dB_arr = [52:0.2:53];
rho_max_arr = 10.^(rho_max_dB_arr/10);
rho_s_arr = rho_max_arr - rho_c; 

OPT_alpha_mat = zeros(N,length(rho_max_dB_arr)); % N x number of rho

for i=1:length(rho_max_dB_arr)
    rho_s = rho_s_arr(i);
    % Dinkelbach’s Method for the fractional form and Successive Convex Approximation (SCA) for the non-convex subproblems
    d_vec = beta_AR; % AP to Repeater channel
    v = sqrt(beta_ADR).*exp(-1i*2*pi*fc*tau_ADn'); % Repeater to Target channel
    % Constants
    c1 = sigma_rcs2 * beta_AD * (rho_c + rho_s * M);
    C2 = sigma_rcs2 * (rho_c + rho_s * M);
    % 2. Optimization Setup
    max_iter = 20;          % Max Dinkelbach iterations
    tol = 1e-8;             % Convergence tolerance
    alpha_opt = alpha_max.*ones(N,1)./2; %* rand(N,1); % Initial feasible point
    eta = 0;                % Initial SINR value (Dinkelbach parameter)
    % Pre-calculate constant matrices for speed
    % Numerator term: d.' * B * v = q.'*alpha + alpha.'*W*alpha 
    q = d_vec .* v;
    W = diag(d_vec) * H_RR * diag(v);

    alpha_opt = compute_alpha_opt(d_vec, v, H_RR, alpha_max, ...
                                       sigma_rcs2, rho_c, rho_s, M, ...
                                       sigma_R2, sigma_AP2);
    OPT_alpha_mat(:,i) = alpha_opt;
end
% repeater selection
repeater_masking = abs(OPT_alpha_mat) > alpha_max*1e-2;
OPT_alpha_mat_masked = OPT_alpha_mat.*repeater_masking*N;

%% Evaluate empirical Sensing SINR over N repeaters via solution ALPHAs above

mc=5000;  %temp1=0; temp2=0; temp3 = 0; temp4=0; temp5 =0;
rho_max_dB_arr = [52:0.2:53];
rho_max_arr = 10.^(rho_max_dB_arr/10);
rho_s_arr = rho_max_arr - rho_c;


sensing_SINR_arr = zeros(1,length(rho_max_arr));
sensing_SINR_arr_ins = zeros(1,length(rho_max_arr));

for j=1:length(rho_max_arr)
    gamma_sen_sim_nu = 0;  gamma_sen_sim_denom = 0; SINR_ins = zeros(mc,1);
    
    % selection ALPHAs
    ALPHA_util = diag(OPT_alpha_mat(:,j)); % from previous section
    %ALPHA_util = eye(N).*alpha_max; % set ALPHA as MAX alpha
    %ALPHA_util = eye(N).*alpha_max*0; % set ALPHA=0 without repeaters 
    % masking
    %ALPHA_util = diag(OPT_alpha_mat_masked(:,j)); % masked alpha x N
    
    rho_s = rho_s_arr(j);
    for i=1:mc
        h_au = sqrt(beta_AU/2).*(randn(M,1)+1i*randn(M,1));  % randomize user channel form channel estimate
        A = conj(h_au)./norm(h_au); B = conj(a(0))./norm(a(0));
        U = [A,B];
        w_s = (eye(M) - U/(U'*U)*U')*conj(w_s_p) ;
        % w_s = conj(w_s_p); % perform veary bad in sensing SINR 
        w_s = (w_s)./norm(w_s); % normalized sensing precoder
        A0_n = conj(a(0))*conj(a(0))'./norm(a(0))^2 ;
        w_c = (eye(M) - A0_n)*conj(h_au);
        % w_c = conj(h_au); % perform veary bad in sensing SINR 
        w_c = (w_c)./norm(w_c); % normalzied communication precoder
        s_s = qammod(randi([0,3],1),4,'UnitAveragePower',true); s_c = qammod(randi([0,3],1),4,'UnitAveragePower',true); %s_s=1; s_c=1;
        x = sqrt(rho_s)*w_s*s_s + sqrt(rho_c)*w_c*s_c ;
        noise_R = sqrt(sigma_R2/2)*(randn(N,1)+1i*randn(N,1));
        noise_AP = sqrt(sigma_AP2/2)*(randn(M,1)+1i*randn(M,1));
        noise_UE = sqrt(sigma_AP2/2)*(randn()+1i*randn());
        rcs_reli = sqrt(sigma_rcs2);%*(randn(1,1)+1i*randn(1,1));

        numer_ins = norm(rcs_reli.*H_AD*x + H_AR/(eye(N) + ALPHA_util*H_RR)*ALPHA_util*conj(H_ADR')*x.*rcs_reli)^2;
        gamma_sen_sim_nu = gamma_sen_sim_nu + numer_ins; % numerator without repeater interaction
        % gamma_sen_sim_nu = gamma_sen_sim_nu + norm(H_AR*(eye(N) + ALPHA*H_RR)*ALPHA*conj(H_ADR')*x)^2; % numerator without repeater interaction
        
        % gamma_sen_sim_denom = gamma_sen_sim_denom + norm(H_AR*ALPHA*(conj(H_AR)'*x + noise_R) + noise_AP)^2 ; 
        denom_ins = norm(H_AR/(eye(N) + ALPHA_util*H_RR)*ALPHA_util*(conj(H_AR')*x + noise_R) + noise_AP)^2;
        gamma_sen_sim_denom = gamma_sen_sim_denom + denom_ins ; 
        %gamma_sen_sim_denom = gamma_sen_sim_denom + norm(H_AR*(eye(N) + ALPHA*H_RR)*ALPHA*( noise_R) + noise_AP)^2 ; 
        SINR_ins(i) = numer_ins/denom_ins;
        %temp1 = temp1 + norm(H_AD*x)^2; 
        % noise_only = noise_only + norm(noise_AP)^2; 
        %temp2 = temp2 + norm((H_AR*ALPHA)*conj(H_ADR)'*x)^2; 
        %temp3 = temp3 + norm((H_AR/(eye(N)-ALPHA*H_RR)*ALPHA)*conj(H_ADR)'*x)^2;
        %temp4 = temp4 + abs(conj(h_au')*sqrt(rho_s)*w_s)^2;
        %temp5 = temp5 + abs(rcs_reli)^2;
    end
    sensing_SINR_arr(j) = gamma_sen_sim_nu/gamma_sen_sim_denom;
    sensing_SINR_arr_ins(j) = mean(SINR_ins); 
    j
end


%% 
%figure 
semilogy(rho_max_dB_arr,sensing_SINR_arr,'--r',LineWidth=2) 
%legend('${\alpha}=\alpha_{\textrm{opt}}, N=5$','$\alpha=\alpha_{\textrm{max}}, N=5$','without repeaters','Interpreter','latex','FontSize', 14,'Location', 'southeast') 
legend('${\alpha}=\alpha_{\textrm{opt}}, N=5$','$\alpha=\alpha_{\textrm{max}}, N=5$','without repeaters','${\alpha}=\alpha_{\textrm{opt}}, N=10$','Interpreter','latex','FontSize', 14,'Location', 'southeast') 

hold on 
grid on 
xlabel('$\rho_{max}$ [dBm]','Interpreter','latex','FontSize', 18),ylabel('$\gamma_s$','Interpreter','latex','FontSize', 18)
%% 
%figure 
semilogy(rho_max_dB_arr,sensing_SINR_arr,':b',LineWidth=2) 
%legend('${\alpha}=\alpha_{\textrm{opt}}, N=5$','$\alpha=\alpha_{\textrm{max}}, N=5$','without repeaters','Interpreter','latex','FontSize', 14,'Location', 'southeast') 
legend('${\alpha}=\alpha_{\textrm{opt}}, N=10$','$\alpha_{\hat{n}}=N\alpha_{\textrm{max}}, \textrm{1-active}$','without repeaters','Interpreter','latex','FontSize', 14,'Location', 'southeast') 

hold on 
grid on 
xlabel('$\rho_{\textrm{max}}$ [dBm]','Interpreter','latex','FontSize', 18),ylabel('$\gamma_s$','Interpreter','latex','FontSize', 18)



%%  plot position
figure
scatter(AP_pos(1,:),AP_pos(2,:),'r')
hold on
scatter(D_pos(1,:),D_pos(2,:),'k')
scatter(Repeater_pos(1,:),Repeater_pos(2,:),'b')

%AP_pos = [0;0]; D_pos = [l_AD*cos(theta_drone) ; l_AD*sin(theta_drone)]; Repeater_pos
