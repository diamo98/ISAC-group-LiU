function [H_all,H_all_fullloss,H01,H02,H03] = Gen_Hi(M,num_points_tot,h_pi_2_aprx,h_ap1_2_pi,h_ap2_2_pi,h_ap3_2_pi,beta_j_Pq,beta_j_Pq_fullloss,radius_tx2points,radius_points2rx,fc, h_direct_dep, h_direct_ari,AP_positions,lambda)

    H1 = zeros(M,M,num_points_tot);
    H2 = zeros(M,M,num_points_tot);
    H3 = zeros(M,M,num_points_tot);
    bistatic_range = radius_tx2points+radius_points2rx;
    light_speed = 3e8;

    for i=1:num_points_tot 
        %H1(:,:,i) = h_pi_2_aprx(:,i)*conj(h_ap1_2_pi(:,i)');
        %H2(:,:,i) = h_pi_2_aprx(:,i)*conj(h_ap2_2_pi(:,i)');
        %H3(:,:,i) = h_pi_2_aprx(:,i)*conj(h_ap3_2_pi(:,i)');
        % with delay
        H1(:,:,i) = h_pi_2_aprx(:,i)*conj(h_ap1_2_pi(:,i)').*exp(-1i*2*pi*fc*bistatic_range(1,i)/light_speed); % .*exp(-1i*2*pi*fc*t)
        H2(:,:,i) = h_pi_2_aprx(:,i)*conj(h_ap2_2_pi(:,i)').*exp(-1i*2*pi*fc*bistatic_range(2,i)/light_speed);
        H3(:,:,i) = h_pi_2_aprx(:,i)*conj(h_ap3_2_pi(:,i)').*exp(-1i*2*pi*fc*bistatic_range(3,i)/light_speed);
        % H = [H1,H2,H3]; % [Mr x J*M]
    end

    Q = num_points_tot;
    H_all = zeros(M,M*3,Q); %J=3
    H_all_fullloss = zeros(M,M*3,Q); %J=3

    for i=1:Q
        %H1_all(:, ((i-1)*M+1):((i-1)*M+M) ) = H1(:,:,i);
        H_all(:,:,i) = [H1(:,:,i).*sqrt(beta_j_Pq(1,i)), H2(:,:,i).*sqrt(beta_j_Pq(2,i)), H3(:,:,i).*sqrt(beta_j_Pq(3,i))];
        H_all_fullloss(:,:,i) = [H1(:,:,i).*sqrt(beta_j_Pq_fullloss(1,i)), H2(:,:,i).*sqrt(beta_j_Pq_fullloss(2,(i))), H3(:,:,i).*sqrt(beta_j_Pq_fullloss(3,i))];
    end
    
    direct_r1 = norm(AP_positions(:,end)-AP_positions(:,1));
    direct_r2 = norm(AP_positions(:,end)-AP_positions(:,2));
    direct_r3 = norm(AP_positions(:,end)-AP_positions(:,3));
    
    %H01 = h_direct_ari(:,1)*conj(h_direct_dep(:,1)').*sqrt(lambda^2/(4*pi*direct_r1)^2).*exp(-1i*2*pi*fc*direct_r1/light_speed);
    %H02 = h_direct_ari(:,2)*conj(h_direct_dep(:,2)').*sqrt(lambda^2/(4*pi*direct_r2)^2).*exp(-1i*2*pi*fc*direct_r2/light_speed);
    %H03 = h_direct_ari(:,3)*conj(h_direct_dep(:,3)').*sqrt(lambda^2/(4*pi*direct_r3)^2).*exp(-1i*2*pi*fc*direct_r3/light_speed);
    H01 = h_direct_ari(:,1)*conj(h_direct_dep(:,1)').*sqrt(1/(direct_r1)^2).*exp(-1i*2*pi*fc*direct_r1/light_speed);
    H02 = h_direct_ari(:,2)*conj(h_direct_dep(:,2)').*sqrt(1/(direct_r2)^2).*exp(-1i*2*pi*fc*direct_r2/light_speed);
    H03 = h_direct_ari(:,3)*conj(h_direct_dep(:,3)').*sqrt(1/(direct_r3)^2).*exp(-1i*2*pi*fc*direct_r3/light_speed);

end