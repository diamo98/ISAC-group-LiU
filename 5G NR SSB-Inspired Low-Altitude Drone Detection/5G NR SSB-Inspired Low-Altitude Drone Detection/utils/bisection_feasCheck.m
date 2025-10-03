function [X_all,rho_ssb_all] = bisection_feasCheck(len_mat,J,M,slack_t_db,ebsilon,H_all,V,precoderSSB_GoB,rcs_var, ... 
    beta_clutter,beta_j_ue,noise_var,gamma_ue,tier1_rho,Pmax,h_direct_link)

% store solution X, rho, slack 
slack_t_bisection = zeros(1,len_mat); X_all = zeros(J*M,J*M,len_mat); rho_ssb_all = zeros(J,len_mat);
block = slack_t_db;  % bisection params

for i=1:len_mat % iterate each point q
    
A = H_all(:,:,i)'*V(:,i)*V(:,i)'*H_all(:,:,i); % JM x JM 
f = precoderSSB_GoB(:,i); % not yet conjugate 
zeros_M = zeros(M,M); n=J*M; iteration_cnt = 1;

while block > ebsilon || (~isFeasible) % bisection loop
    
    [isFeasible, solution_X,solution_rho] = checkFeasibilityCVX(A,rcs_var, 10^(slack_t_db/10),beta_clutter,zeros_M,M, ... 
        Pmax,f,beta_j_ue,noise_var,gamma_ue,n,J,tier1_rho,slack_t_db, h_direct_link);
    
    while iteration_cnt == 1 && (isFeasible) % make sure at the starting slack t should be infeasible
        slack_t_db = slack_t_db*2;
        block = slack_t_db;
        [isFeasible, solution_X,solution_rho] = checkFeasibilityCVX(A,rcs_var, 10^(slack_t_db/10),beta_clutter,zeros_M,M,Pmax,f,beta_j_ue,noise_var,gamma_ue,n,J,tier1_rho,slack_t_db);
    end

    if (isFeasible) % feasible
        slack_t_db_pre = slack_t_db; % store feasible slack in the final case
        slack_t_db = slack_t_db + block/2 % print slack
        block = block/2;
    else % infeasible
        slack_t_db = slack_t_db - block/2 % print slack
        block = block/2;
    end
    iteration_cnt = iteration_cnt+1;
end % while bisection
    slack_t_bisection(i) = slack_t_db_pre; % 1 x len_mat (SCIR)
    X_all(:,:,i) = solution_X;             % JM x JM x len_mat
    rho_ssb_all(:,i) = solution_rho;           % J x len_mat
end % i

end