function [precoder_vec,s_pwr_mat, combine_s_pwr_p_ssb] = X_2_precoder_pwr_scale(M,J,len_point,X_all,rho_ssb_all,Pmax)
    precoder_vec =  zeros(M,J,len_point); s_pwr_mat = zeros(J,len_point); 
    for i=1:len_point
    [vec_i, val_i] = eig(X_all(:,:,i)); % eigen vec and val 
    eigenvalues = diag(val_i); % Extract diagonal of D as a vector
    [sorted_eigenvalues, idx] = sort(eigenvalues); % Sort eigenvalues in ascending order
    sorted_eigenvectors = vec_i(:, idx); % Reorder eigenvectors accordingly
    precoder_Pi = sorted_eigenvectors(:,end); % eig vec corresponds to dom eig value
    %s_pwr_point_i = [trace([eye(M),zeros_M,zeros_M]*X_all(:,:,i)*[eye(M),zeros_M,zeros_M]'), trace([zeros_M,eye(M),zeros_M]*X_all(:,:,i)*[zeros_M,eye(M),zeros_M]'), trace([zeros_M,zeros_M,eye(M)]*X_all(:,:,i)*[zeros_M,zeros_M,eye(M)]') ];
    precoder_Pi_mat = [precoder_Pi(1:M),precoder_Pi(M+1:2*M),precoder_Pi(2*M+1:3*M)];
    % sensing pwr scaling with optimizer
    %precoder_Pi_pwr_scale = precoder_Pi_mat./([norm(precoder_Pi_mat(:,1)), norm(precoder_Pi_mat(:,2)), norm(precoder_Pi_mat(:,3))]) .*sqrt(s_pwr_point_i);
    % sensing pwr scaling with total pwr budget
    precoder_Pi_pwr_scale = precoder_Pi_mat./([norm(precoder_Pi_mat(:,1)), norm(precoder_Pi_mat(:,2)), norm(precoder_Pi_mat(:,3))]) .*sqrt( Pmax-rho_ssb_all(:,i)' );
    
    precoder_vec(:,:,i) = precoder_Pi_pwr_scale;
    s_pwr_mat(:,i) = [norm(precoder_vec(:,1,i))^2; norm(precoder_vec(:,2,i))^2; norm(precoder_vec(:,3,i))^2];
    end
    combine_s_pwr_p_ssb = s_pwr_mat + rho_ssb_all; % p tx at each AP at each time slot
end

