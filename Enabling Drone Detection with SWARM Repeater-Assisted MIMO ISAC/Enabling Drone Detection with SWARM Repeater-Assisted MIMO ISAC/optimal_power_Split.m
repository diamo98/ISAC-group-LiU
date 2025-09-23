function [rho_c, rho_s] = optimal_power_Split(gamma_UE,M,sigma_ue2,beta_AU,rho_max)
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
end
