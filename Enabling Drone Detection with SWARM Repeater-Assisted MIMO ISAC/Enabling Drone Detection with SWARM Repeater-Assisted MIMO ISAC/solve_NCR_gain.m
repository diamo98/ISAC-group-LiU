function [alpha, rho_c, rho_s, gamma_opt, info] = solve_NCR_gain( ...
    beta_AD, beta_AU, beta_A, beta_AD_n, ...
    sigma_r2, sigma_ap2, sigma_ue2, gamma_UE, ...
    M, rho_max, alpha_max, opts)
% SOLVE_PROBLEM  Closed-form power split + Dinkelbach for alphas (global optimum).
% Inputs follow the statement; vectors are Nx1. opts.tol (1e-12), opts.maxit (1000), opts.verbose (0/1).
% Outputs:
%   alpha[Nx1]  : optimal {alpha_n} in {0, alpha_max}
%   rho_c,rho_s : closed-form optimal powers
%   gamma_opt   : optimal objective value
%   info        : diagnostics

    if nargin < 12 || isempty(opts), opts = struct; end
    tol     = getOpt(opts,'tol',1e-15);
    maxit   = getOpt(opts,'maxit',1000);
    verbose = getOpt(opts,'verbose',0);

    beta_A    = beta_A(:);
    beta_AD_n = beta_AD_n(:);
    N = numel(beta_A);
    assert(numel(beta_AD_n)==N, 'beta_A and beta_AD_n must match.');
    assert(M > 0 && beta_AU > 0, 'M and beta_AU must be positive.');
    assert(sigma_ap2 >= 0 && sigma_r2 >= 0 && sigma_ue2 >= 0, 'noise powers >= 0');

    % ==== 1) Closed-form optimal (rho_c, rho_s) ====
    a = gamma_UE / M;
    b = a * (sigma_ue2 / beta_AU);
    if rho_max < b - 1e-14
        rho_c = rho_max;
        rho_s = 0;
        %alpha = zeros(N,1);
        %rho_c = NaN; rho_s = NaN; gamma_opt = -Inf;
        %info = struct('feasible',false,'reason', ...
        %    'Infeasible: rho_max < (gamma_UE/M)*(sigma_ue^2/beta_AU).');
        %if verbose, warning(info.reason); end
        %return;
    else
        rho_s = max(0, (rho_max - b) / (1 + a));
        rho_c = rho_max - rho_s;
    end
    Lstar = rho_c + M * rho_s;

    % ==== 2) Dinkelbach for alpha-subproblem ====
    A = beta_AD;     B = sigma_ap2;
    p = beta_A .* beta_AD_n;       % numerator coeffs
    r = sigma_r2 .* beta_A;        % denominator coeffs
    u = alpha_max^2;

    if u <= 0
        alpha = zeros(N,1);
        z = A / max(B, eps);
        gamma_opt = Lstar * z;
        info = packInfo(true, 0, z, [], rho_c, rho_s, Lstar);
        return;
    end

    z = max(A / max(B, eps), 0);   % safe init
    it = 0; fval = inf; mask = false(N,1);

    while it < maxit
        it = it + 1;

        % on/off choice at this z
        s = p - z .* r;
        mask = (s > 0);
        t = u * double(mask);      % t_n = alpha_n^2

        % update z and residual
        num = A + sum(p .* t);
        den = B + sum(r .* t); den = max(den, eps);
        fval = num - z * den;
        zNew = num / den;

        if verbose
            fprintf('It=%3d, z=%.12g, res=%.3e, |dz|=%.3e, active=%d\n', ...
                it, z, fval, abs(zNew - z), nnz(mask));
        end

        if abs(fval) <= tol || abs(zNew - z) <= tol * max(1,abs(z))
            z = zNew; break;
        end
        z = zNew;
    end

    alpha = sqrt(t);
    gamma_opt = Lstar * z;
    info = packInfo(true, it, z, find(mask), rho_c, rho_s, Lstar);
end

% ---- helpers ----
function v = getOpt(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
function info = packInfo(feasible, it, z, active_idx, rho_c, rho_s, Lstar)
    info = struct('feasible',feasible,'iters',it,'z_ratio',z, ...
                  'active_idx',active_idx,'rho_c',rho_c,'rho_s',rho_s,'Lstar',Lstar);
end
