clear; clc; close all;


%%
c = physconst("LightSpeed");
options = optimoptions('fmincon','display','none',... 
    'MaxFunctionEvaluations', 1e4, 'EnableFeasibilityMode', true);

MC_RUNS = 3;

parameters = [6,7];
n_params = numel(parameters);

g_mc = cell(MC_RUNS, n_params);
crb_mc = cell(MC_RUNS, n_params);
P_mc = cell(MC_RUNS, n_params);
convergence_mc = cell(MC_RUNS, n_params);
ue_sinr_mc = cell(MC_RUNS, n_params);
sigma_mc = cell(MC_RUNS, n_params);
alpha_mc = cell(MC_RUNS, n_params);
w_mc = cell(MC_RUNS, n_params);

for param_idx = 1:n_params
    tmp = parameters(param_idx);
    for mc = 1:MC_RUNS
        tic
        tmp
        mc
        convergence = 1;

        % Distances
        user_repeater_d = 20;
        d = 100;
        repeater_AP_d = d;
        
        % Path loss exponents
        fspl_exp = 2;
        pl_exp = 3;    
        
        % Variances
        noise_var = 10^(-94/100)/1000;
        channel_est_var = 10^(-20/10);
        
        sH = channel_est_var;
        sz = noise_var * (1/(repeater_AP_d^pl_exp));
        se = noise_var;
        sh = (1/(user_repeater_d^pl_exp));
        sg = (1/(user_repeater_d^pl_exp)) * (1/(repeater_AP_d^pl_exp));
        
        % RCS
        rcs_var = 10^(10/10);
        sigma = sqrt(rcs_var/2)*(randn(1,1) + 1j*randn(1,1));
        
        % Misc
        M = 64;
        N = 128;
        Ns = 128;
        delta_f = 120*10^3;
        phi = 30*(pi/180);
        a_phi = exp((1j*(0:M-1)*pi*cos(phi))).';
        g_k = sqrt(sg/2)*(randn(M,Ns) + 1j*randn(M,Ns));
        
        % Constraints
        P_max = M * (10^(35/10)/1000);
        UE_SINR_min = 10^(-2/10);
        
        % Thresholds
        alpha_threshold = 1;
        w_threshold = P_max/M;
        max_its = 20;
        
        % Initial values
        alpha = 1;
        w = (1/sqrt(2) + 1j/sqrt(2)) * ones(M,1);
        
        it = 1;
        threshold = true;
        
        while threshold
            alpha_prev = alpha;
            w_prev = w;
            %it
        
        
            fun = @(x) objfun(x, d, sigma, M, N, a_phi, sH, sz, se, Ns, delta_f, c);
        
            nonlcon = @(x) nonlconfun(x, g_k, Ns, se, sh, UE_SINR_min, P_max, M);
            x0 = [real(w); imag(w); alpha];
            x = fmincon(fun,x0,[],[],[],[],[],[],nonlcon, options);
            w = (x(1:M) + 1j*x(M+1:end-1));
            alpha = x(end);
        
            crb = objfun(x, d, sigma, M, N, a_phi, sH, sz, se, Ns, delta_f, c);
        
            % convergence checks
            threshold = ...
                (norm(w-w_prev) > w_threshold) ...
                || (abs(alpha - alpha_prev) > alpha_threshold);
        
            % if not successful, re-randomize and re-run (is this stupid?)
            constraints = nonlconfun(x, g_k, Ns, se, sh, UE_SINR_min, P_max, M);
            if (~threshold) && (~all(constraints < 0) || crb < 0)
                threshold = true;
                sigma = sqrt(rcs_var/2)*(randn(1,1) + 1j*randn(1,1));
                w = (1/sqrt(2) + 1j/sqrt(2)) * ones(M,1);
                g_k = sqrt(sg/2)*(randn(M,Ns) + 1j*randn(M,Ns)); % each column is one k
                alpha = 1;
            end
            if it >= max_its
                if ~all(constraints < 0)
                    convergence = 0;
                end
                break;
            end
            it = it + 1;
        end
        toc
        g_mc{mc, param_idx} = g_k;
        crb_mc{mc, param_idx} = crb;
        P_mc{mc, param_idx} = norm(w)^2;
        convergence_mc{mc, param_idx} = convergence;
        ue_sinr_mc{mc, param_idx} = (alpha^2 * sum(abs(g_k.'*w).^2))/(Ns*x(end)^2 * sh*se + Ns*se);
        sigma_mc{mc, param_idx} = sigma;
        alpha_mc{mc, param_idx} = alpha;
        w_mc{mc, param_idx} = w;
    end
end


clear fun nonlcon
save("test_run.mat");
exit

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
function val = objfun(x, d, sigma, M, N, a_phi, sH, sz, se, Ns, delta_f, c)
    alpha = x(end);
    w = (x(1:M) + 1j*x(M+1:end-1));

    sigma_R = real(sigma);
    sigma_I = imag(sigma);
    psi = (M*N*real(w'*(a_phi*a_phi')*w))/(alpha^2 * sH*real(w'*w) + alpha^2 * sz + se);
    C = (Ns/(d^2)) + (16*pi^2 * delta_f^2 * Ns*(Ns-1)*(2*Ns-1))/(6*c^2);
    S_R = (-sigma_R*Ns)/d - (sigma_I*4*pi*delta_f*Ns*(Ns-1))/(2*c);
    S_I = (-sigma_I*Ns)/d + (sigma_R*4*pi*delta_f*Ns*(Ns-1))/(2*c);
    
    val = 8*d^6 / (psi^3*(abs(sigma)^2 *C - S_R - S_I));
    %val = d^12 / (psi^3*(abs(sigma)^2 *C - S_R - S_I));
end