function [X_all,rho_ssb_all] = bisection_feasCheck_speed(len_point,J,M,slack_t_db,ebsilon,H_all,V,precoderSSB_GoB,rcs_var, ... 
    beta_clutter,beta_j_ue,noise_var,gamma_ue,Pmax,h_direct_dep)

% NOTE that this function is designed for only the case with hexagonal arrangment.

X_all = zeros(J*M,J*M,len_point); rho_ssb_all = zeros(J,len_point);
block = slack_t_db;  % bisection params

for i=1:len_point % iterate each point q
    
    A = H_all(:,:,i)'*V(:,i)*V(:,i)'*H_all(:,:,i); % JM x JM 
    f = precoderSSB_GoB(:,i); % not yet conjugate 
    zeros_M = zeros(M,M); JM=J*M; iteration_cnt = 1; isFeasible = false;
    
    while (~isFeasible) % bisection loop
        
        [isFeasible, solution_X,solution_rho] = checkFeasibilityCVX(A,rcs_var, 10^(slack_t_db/10),beta_clutter, ... 
            zeros_M,M,Pmax,f,beta_j_ue,noise_var,gamma_ue,JM,J,slack_t_db,h_direct_dep);
        
    
        if (~isFeasible) % infeasible
            slack_t_db = slack_t_db - block/2; % print slack 
            %block = block/2; 
            slack_t_db
        end
        iteration_cnt = iteration_cnt+1;
    end % while bisection
        %slack_t_bisection(i) = slack_t_db_pre; % 1 x len_mat (SCIR)
        X_all(:,:,i) = solution_X;             % JM x JM x len_mat
        rho_ssb_all(:,i) = solution_rho;           % J x len_mat
        fprintf('Run i=%d out of %d \n',i,len_point)
end % i

end