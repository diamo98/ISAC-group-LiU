function alpha_opt = compute_alpha_opt(d, v, H_RR, alpha_max, ...
                                       sigma_rcs2, rho_c, rho_s, M, ...
                                       sigma_R2, sigma_AP2)

    % Number of repeaters
    N = length(d);

    % -------------------------------------------------------------
    % Build matrices Q and R such that:
    % numerator   = c1 + x' * Q * x
    % denominator = sigma_R^2 * x' * R * x + sigma_AP^2
    % where x = alpha .* d
    % -------------------------------------------------------------

    % c1 constant term
    c1 = sigma_rcs2 * (rho_c + rho_s*M);

    % Build Q and R explicitly
    % Q corresponds to sensing gain:
    %   d^T * B * v
    % with B = (I + Phi_a H_RR) Phi_a
    %
    % After reformulation: numerator(x) = x.' * Q * x

    % Build effective matrices
    D = diag(d);

    % Construct H1 and H2 mapping x = Phi_a * d
    H1 = eye(N);                   % from Phi_a
    H2 = H_RR;                     % from Phi_a * H_RR * Phi_a

    % Define the "v" term after diagonal multiplication
    Vd = diag(v);

    % Q is built from:
    %   g^T v + g^T H_RR (v ⊙ alpha)
    % in quadratic form on x = alpha ⊙ d

    Q = sigma_rcs2 * (rho_c + rho_s*M) * ...
        ( D*Vd*D + D*H_RR*Vd*D + D*Vd*H_RR*D + D*H_RR*Vd*H_RR*D );

    % R is built from:
    %   || g + Phi_a H_RR g ||^2
    % = (g + H_RR g)' * (g + H_RR g)
    %
    % In quadratic form: x' R x

    R = (D + D*H_RR)' * (D + D*H_RR);

    % -------------------------------------------------------------
    % Solve generalized eigenvalue problem:
    %   A x = lambda C x
    % where A = Q, C = sigma_R2 * R
    % -------------------------------------------------------------

    A = Q;
    C = sigma_R2 * R;

    % Compute principal generalized eigenvector
    [vec, val] = eigs(A, C, 1, 'largestreal');
    x_star = vec;   % optimal x

    % -------------------------------------------------------------
    % Map back to alpha = x ./ d
    % -------------------------------------------------------------
    alpha_unscaled = x_star ./ d;

    % -------------------------------------------------------------
    % Apply |alpha_n| <= alpha_max
    % -------------------------------------------------------------
    scale = min(alpha_max ./ abs(alpha_unscaled));
    alpha_opt = scale * alpha_unscaled;

end