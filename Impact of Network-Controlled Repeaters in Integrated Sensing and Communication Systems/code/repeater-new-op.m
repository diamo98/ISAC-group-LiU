clear; clc; close all;


%%
c = physconst("LightSpeed");
options = optimoptions('fmincon','display','none',... 
    'MaxFunctionEvaluations', 1e4, 'EnableFeasibilityMode', true);...
    %'PlotFcn','optimplot');

%123
% Distances
user_repeater_d = 100;

% using same distance to target and repeater
d = 500;
repeater_AP_d = d;

fspl_exp = 2;
pl_exp = 4;

% Variances
noise_var = 10^(-94/10)/1000; % -94 dBm
channel_est_var = 10^(-20/10); %anywhere from 0 to -40 dB variance?
%https://www.researchgate.net/figure/Simulation-for-the-variance-of-the-channel-estimation-error-in-Gaussian-noise_fig5_292314336

sH = channel_est_var;
sz = noise_var * (1/(repeater_AP_d^pl_exp));
se = noise_var;
sh = (1/(user_repeater_d^pl_exp));
sg = (1/(user_repeater_d^pl_exp)) * (1/(repeater_AP_d^pl_exp));

% RCS
rcs_var = 10^(5/10); % 5 dB
sigma = sqrt(rcs_var/2)*(randn(1,1) + 1j*randn(1,1));
%sigma = sqrt(rcs_var/2) + 1j*sqrt(rcs_var/2);

% Misc
M = 64;
N = 256;
Ns = 256;
delta_f = 120*10^3; % 120 kHz
phi = 30*(pi/180);
phi2 = 75*(pi/180);
a_phi = exp((1j*(0:M-1)*pi*cos(phi))).';
a_phi2 = exp((1j*(0:M-1)*pi*cos(phi2))).';
g_k = sqrt(sg/2)*(randn(M,Ns) + 1j*randn(M,Ns)); % each column is one k
%g_k = sqrt(sg/2)*a_phi2.*ones(M,Ns);

% Constraints
P_max = (10^(45/10)/1000); % 35 dBm gNodeB har 43-55 dbm https://www.ericsson.com/en/blog/2023/8/breaking-the-energy-curve
% g8 kanske med 45-55? fix alpha kan behöva vara 20dB då?
UE_SINR_min = 10^(2/10); % 3 dB

% Thresholds
alpha_threshold = 1;
w_threshold = P_max/M;
max_its = 20;

% Initial values
alpha = 1;%10^(40/10); % 30dB
w = (1/sqrt(2) + 1j/sqrt(2)) * ones(M,1);
%w = sqrt(P_max)*(steerings./(norm(steerings)));


A = @(w) (2*M*N*real(w'*(a_phi*a_phi')*w));
B = @(alpha, w) (alpha^2 * sH*real(w'*w) + alpha^2 * sz + se);


it = 1;
threshold = true;
tic
while threshold
    alpha_prev = alpha;
    w_prev = w;
    %it

    %psi = (2*M*N*real(w'*(a_phi*a_phi')*w))/(alpha^2 * sH*real(w'*w) + alpha^2 * sz + se);

    if mod(it, 2) == 0
        cvx_clear
        cvx_begin
            variable w(M,1) complex
            variable y(1,1)
            obj = 2*y*A(w) - y^2 * B(w)^2;

        maximize obj
        subject to
            numerator = (alpha^2 * sum(abs(g_k.'*w).^2));
            denominator = (Ns*x(end)^2 * sh*se + Ns*se);
            numerator/denominator >= UE_SINR_min;
    else
        %variable tau_n nonnegative
        %variable tau_d nonnegative
    end

    crb = objfun(w, alpha, d, sigma, M, N, a_phi, sH, sz, se, Ns, delta_f, c);

    % convergence checks
    % threshold = ...
    %     (norm(w-w_prev) > w_threshold) ...
    %     || (abs(alpha - alpha_prev) > alpha_threshold);

    if it >= max_its
        % if constraints... etc
        disp("did not converge")
        break;
    end
    it = it + 1;
end
toc
disp(sprintf("std crb: %2.9f", sqrt(crb)));
disp(sprintf("power: %2.0f", norm(w)^2));
disp(sprintf("alpha db: %2.0f", 10*log10(alpha)));
disp(sprintf("ue sinr: %2.2f", (alpha^2 * sum(abs(g_k.'*w).^2))/(Ns*alpha^2 * sh*se + Ns*se)));
%%
figure
angles = linspace(0,pi/2,1000);
beams = [];
for a = angles
    beams = [beams, abs(w'*exp((1j*(0:M-1)*pi*cos(a))).'*exp((1j*(0:M-1)*pi*cos(a)))*w)];
end
plot(angles*180/pi, beams);
title("w");

figure
hold on
plot(real(w), imag(w), 'ro');
hold off



%% 
% non linear constraints
function [c,ceq] = nonlconfun(x, g_k, Ns, se, sh, UE_SINR_min, P_max, M)
    alpha = x(end);
    w = (x(1:M) + 1j*x(M+1:end-1));

    c(1) = UE_SINR_min - (alpha^2 * sum(abs(g_k.'*w).^2))...
        /(Ns*x(end)^2 * sh*se + Ns*se);
    c(2) = norm(w)^2 - P_max; % tri ineq?
    %c(2) = sum(abs(w(1:M) + w(M+1:end)).^2) - P_max;
    c(3) = -alpha;
    ceq = [];
end

% objective function
function val = objfun(w, alpha, d, sigma, M, N, a_phi, sH, sz, se, Ns, delta_f, c)
    sigma_R = real(sigma);
    sigma_I = imag(sigma);
    psi = (2*M*N*real(w'*(a_phi*a_phi')*w))/(alpha^2 * sH*real(w'*w) + alpha^2 * sz + se);
    C = (4*Ns/(d^2)) + (16*pi^2 * delta_f^2 * Ns*(Ns-1)*(2*Ns-1))/(6*c^2);
    S_R = (-sigma_R*2*Ns)/d - (sigma_I*4*pi*delta_f*Ns*(Ns-1))/(2*c);
    S_I = (-sigma_I*2*Ns)/d + (sigma_R*4*pi*delta_f*Ns*(Ns-1))/(2*c);

    FIM = (psi/d^4).*[abs(sigma)^2 *C, S_R, S_I
                      S_R, Ns, 0
                      S_I, 0, Ns];
    
    val = d^4 / (Ns*psi*(abs(sigma)^2 *Ns*C - S_R^2 - S_I^2));
    %val = d^12 / (psi^3*(abs(sigma)^2 *C - S_R - S_I));
end