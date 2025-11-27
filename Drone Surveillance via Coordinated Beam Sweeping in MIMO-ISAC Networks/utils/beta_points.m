function [beta_j_Pq,beta_j_Pq_fullloss] = beta_points(J,Q,point_positions,AP_positions,lambda) %ignore the constants lambda^2/(4*pi)^3

    [beta_j_Pq] = zeros(J,Q);
    [beta_j_Pq_fullloss] = zeros(J,Q);
    for i=1:Q
        % AP1
        la_1 = norm(point_positions(:,i)-AP_positions(:,1));
        lb = norm(point_positions(:,i)-AP_positions(:,4));
        beta_j_Pq(1,i) = 1/((la_1^2)*(lb^2));
        %beta_j_Pq_fullloss(1,i) = lambda^2/((4*pi)^3*(la_1^2)*(lb^2));
        beta_j_Pq_fullloss(1,i) = beta_j_Pq(1,i);
        % AP2
        la_2 = norm(point_positions(:,i)-AP_positions(:,2));
        beta_j_Pq(2,i) = 1/((la_2^2)*(lb^2));
        %beta_j_Pq_fullloss(2,i) = lambda^2/((4*pi)^3*(la_2^2)*(lb^2));
        beta_j_Pq_fullloss(2,i) = beta_j_Pq(2,i);
        % AP3
        la_3 = norm(point_positions(:,i)-AP_positions(:,3));
        beta_j_Pq(3,i) = 1/((la_3^2)*(lb^2));
        %beta_j_Pq_fullloss(3,i) = lambda^2/((4*pi)^3*(la_3^2)*(lb^2));
        beta_j_Pq_fullloss(3,i) = beta_j_Pq(3,i) ;
        
    end
        
    
    
end
