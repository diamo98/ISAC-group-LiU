function gamma = sensing_gamma(alpha, rho_c, rho_s, M, ...
                               beta_AD, beta_A, beta_AD_n, ...
                               sigma_r2, sigma_ap2)
%SENSING_GAMMA Compute
% gamma = ((rho_s*M + rho_c) * ( beta_AD + sum(alpha_n^2 * beta_A,n * beta_AD,n) )) ...
%          / ( sigma_r2 * sum(alpha_n^2 * beta_A,n) + sigma_ap2 )

    % ensure column vectors
    alpha     = alpha(:);
    beta_A    = beta_A(:);
    beta_AD_n = beta_AD_n(:);

    % basic size check
%    assert(numel(alpha)==numel(beta_A) && numel(beta_A)==numel(beta_AD_n), ...
 %       'alpha, beta_A, and beta_AD_n must have the same length.');

    t   = alpha.^2;
    num = (rho_s*M + rho_c) * ( beta_AD + sum( t .* beta_A .* beta_AD_n ) );
    den = sigma_r2 * sum( t .* beta_A ) + sigma_ap2;

    gamma = num / den;
end