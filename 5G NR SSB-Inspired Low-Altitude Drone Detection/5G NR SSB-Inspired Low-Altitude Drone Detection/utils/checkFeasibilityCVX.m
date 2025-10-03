function [isFeasible, solution_X,solution_rho] = checkFeasibilityCVX(A,rcs_var,slack_t,beta_clutter,zeros_M,M, ... 
    Pmax,f,beta_j_ue,noise_var,gamma_ue,JM,J,tier1_rho,slack_t_db,h_direct_dep)

    % Initialize output
    isFeasible = false;
    solution_X = 0;
    solution_rho = 0;
    
    % Check feasibility using CVX
    %try
    cvx_begin sdp quiet
        variable X(JM,JM) hermitian semidefinite       % Hermitian matrix variable
        variable rho(J)
        maximize real(trace(A * X))    % to only check feasibility
        subject to
            real(trace(A * X))*rcs_var*1e2 >= slack_t*(noise_var + sum(rho)*beta_clutter )*1e2 % *1e2 for tolerance  % 20b
            trace([eye(M),zeros_M,zeros_M]*X*[eye(M),zeros_M,zeros_M]') + rho(1) <= Pmax % 20c
            trace([zeros_M,eye(M),zeros_M]*X*[zeros_M,eye(M),zeros_M]') + rho(2) <= Pmax % 20c
            trace([zeros_M,zeros_M,eye(M)]*X*[zeros_M,zeros_M,eye(M)]') + rho(3) <= Pmax % 20c
            rho(1) >= 0 
            rho(2) >= 0
            rho(3) >= 0
            conj(f)'*[eye(M),zeros_M,zeros_M]*X*[eye(M),zeros_M,zeros_M]'*conj(f) == 0 % 20d 
            conj(f)'*[zeros_M,eye(M),zeros_M]*X*[zeros_M,eye(M),zeros_M]'*conj(f) == 0 % 20d
            conj(f)'*[zeros_M,zeros_M,eye(M)]*X*[zeros_M,zeros_M,eye(M)]'*conj(f) == 0 % 20d
            % direct suppression
            %conj(h_direct_dep(:,1))'*[eye(M),zeros_M,zeros_M]*X*[eye(M),zeros_M,zeros_M]'*conj(h_direct_dep(:,1)) == 0 % 20f
            %conj(h_direct_dep(:,2))'*[zeros_M,eye(M),zeros_M]*X*[zeros_M,eye(M),zeros_M]'*conj(h_direct_dep(:,2)) == 0 % 20f
            %conj(h_direct_dep(:,3))'*[zeros_M,zeros_M,eye(M)]*X*[zeros_M,zeros_M,eye(M)]'*conj(h_direct_dep(:,3)) == 0 % 20f
            
            % SINR ue constraints
            %rho(1)*beta_j_ue(1,1) >= gamma_ue*(rho(1)*sum(beta_j_ue(1,2:3)) + sum(diag(X).*[beta_j_ue(1,1)*ones(M,1); beta_j_ue(1,2)*ones(M,1); beta_j_ue(1,3)*ones(M,1)]) + noise_var)
            %rho(2)*beta_j_ue(2,1) >= gamma_ue*(rho(2)*sum(beta_j_ue(2,2:3)) + sum(diag(X).*[beta_j_ue(2,1)*ones(M,1); beta_j_ue(2,2)*ones(M,1); beta_j_ue(2,3)*ones(M,1)]) + noise_var)
            %rho(3)*beta_j_ue(3,1) >= gamma_ue*(rho(3)*sum(beta_j_ue(3,2:3)) + sum(diag(X).*[beta_j_ue(3,1)*ones(M,1); beta_j_ue(3,2)*ones(M,1); beta_j_ue(3,3)*ones(M,1)]) + noise_var)
            
            %gamma_ue*((beta_j_ue(1,2)*tier1_rho) + noise_var)*10^8  <=  rho(1)*beta_j_ue(1,1)*10^8 % 20e
            %gamma_ue*((beta_j_ue(1,2)*tier1_rho) + noise_var)*10^8  <=  rho(2)*beta_j_ue(2,2)*10^8 % 20e
            %gamma_ue*((beta_j_ue(1,2)*tier1_rho) + noise_var)*10^8  <=  rho(3)*beta_j_ue(3,3)*10^8 % 20e
            gamma_ue*((beta_j_ue(1,2)*rho(1)) + noise_var)*10^8  <=  rho(1)*beta_j_ue(1,1)*10^8 % 20e
            gamma_ue*((beta_j_ue(1,2)*rho(2)) + noise_var)*10^8  <=  rho(2)*beta_j_ue(2,2)*10^8 % 20e
            gamma_ue*((beta_j_ue(1,2)*rho(3)) + noise_var)*10^8  <=  rho(3)*beta_j_ue(3,3)*10^8 % 20e
    cvx_end
        
        if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
            isFeasible = true;
            solution_X = X;
            solution_rho = rho;
            disp('Problem is feasible. Solution found.');
        else
            disp(['Problem is infeasible. CVX status: ' cvx_status]);
        end
        fprintf('Checking feasibility ... ')
end
