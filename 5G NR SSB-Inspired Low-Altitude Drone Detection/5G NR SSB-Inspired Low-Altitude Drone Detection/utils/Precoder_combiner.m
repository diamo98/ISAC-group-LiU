clear all

%%%%%%% setting %%%%%%%%%%%%%
MV=12; MH=12; M=MV*MH; num_aptx = 3; num_aprx = 1; ISD = 500; J=3; T=5;
x=ISD/2; y=ISD/2; z=ISD/2; d=8; rcs_var = 10^(-10/10); noise_var = 10^(-94/10); Pmax = 10^(43/10); %mW      
% define a space volumn [width x height x depth]
radar_vol_width = 20; radar_vol_height = 15; radar_vol_depth = 25;
point_row=ceil(radar_vol_height/d +1); point_col=ceil(radar_vol_width/d +1); point_deep=ceil(radar_vol_depth/d +1);
num_points_tot = point_row*point_col*point_deep; Q = num_points_tot; mod_order = 4;

% define cartisian of points and APs 
[AP_positions,point_positions] = setup_cart(x,y,z,point_row,point_col,point_deep,d);
% define polar between APs_tx and points
[azi_deg_tx2points, ele_deg_tx2points, radius_tx2points] = polar_APs_points(AP_positions,point_positions);
% define polar between points and APs_rx
[azi_deg_points2rx, ele_deg_points2rx, radius_points2rx] = polar_Points_APrx(AP_positions,point_positions);
 
% define H=[H1,H2,..] at point-1 [MxQ]
[h_ap1_2_pi, h_ap2_2_pi, h_ap3_2_pi, h_pi_2_aprx, h_APs_2_pi] = Channel_points2APs(MV,MH, azi_deg_tx2points, ele_deg_tx2points, azi_deg_points2rx, ele_deg_points2rx);
% H1, H2, H3 where [M x M x Q]
[H1,H2,H3,H_all] = Gen_Hi(M,num_points_tot,h_pi_2_aprx,h_ap1_2_pi,h_ap2_2_pi,h_ap3_2_pi);
% combining vector v
V = 1/sqrt(M).*h_pi_2_aprx;

% W= [W1,W2,...,W_R] SSB beams and 
[azimu_GoB_SSB_deg, polar_GoB_SSB_deg, BeamList, nr_tot_SSBbeams, precoderSSB_GoB] = orthogonalGoB_PrecodeGen(MV); % normalized precoder
R = nr_tot_SSBbeams*4;
len_mat = min(Q,R); % Q is set to be less than R
%precoder_Sens_APs = 1./sqrt(M).*conj(h_APs_2_pi); % [M x Q x 3] % to be optimized
precoder_Sens_APs = ones(M,len_mat,3); % [M x Q x 3] % to be optimized

% Gen Comm, SSB Symbols
SSB_symbols1 = qammod(randi([0,mod_order-1],R,1),mod_order); %[4R x 1]
SSB_symbols2 = qammod(randi([0,mod_order-1],R,1),mod_order); %[4R x 1]
SSB_symbols3 = qammod(randi([0,mod_order-1],R,1),mod_order); %[4R x 1]
Sens_symbols = qammod(randi([0,mod_order-1],Q,1),mod_order); %[Q x 1]

% alocating power coef
rho_com = ones(len_mat,1); %rho_sens = ones(J,1);

% beta_j_pq
[beta_j_Pq, beta_ue] = beta_points(J,Q,point_positions,AP_positions);
% beta_j_ue
max_d = x-10; min_d=10;
beta_j_ue =  1./(4*pi*(min_d+rand(len_mat,1).*(max_d-min_d)).^2) ;
max_d_i = 2*x-10; min_d_i=x+10;
beta_i_ue =  1./(4*pi*(min_d_i+rand(len_mat,2).*(max_d_i-min_d_i)).^2);

% Rx signals Y_APr in [M x min{Q,R}] % raw
Y_APr = zeros(M,len_mat);
y_apr_comb = zeros(1,len_mat);
SINR_apr_nom = zeros(len_mat,1);
SINR_apr_denom =zeros(len_mat,1); 
SINR_s =zeros(len_mat,1); 
for i=1:Q % iterate each point
    noise_apr = sqrt(noise_var/2).*( randn(M,len_mat) + 1i*randn(M,len_mat));
    H_ssb_echoes = 1/sqrt(2).*( randn(M,M*T) + 1i*randn(M,M*T));
    nom = 0;
    denom = 0;
    for j=1:J
        rcs = sqrt(rcs_var/2).*(randn(1)+1i*randn(1));
        %Y_APr(:,i) = Y_APr(:,i) + rcs.*sqrt(beta_j_Pq(j,i)).*H_all(:, ((j-1)*M+1):((j-1)*M+M) )*precoder_Sens_APs(:,i,j).*Sens_symbols(i) +...
            %sqrt(rho_com(j)).*(H_ssb_echoes*[precoderSSB_GoB(:,i);precoderSSB_GoB(:,i);precoderSSB_GoB(:,i);precoderSSB_GoB(:,i);precoderSSB_GoB(:,i)]).*SSB_symbols1(i) + ...
            %noise_apr(:,i);
        nom = nom + rcs.*sqrt(beta_j_Pq(j,i)).*H_all(:, ((j-1)*M+1):((j-1)*M+M) )*precoder_Sens_APs(:,i,j).*Sens_symbols(i);
        denom = denom + sqrt(rho_com(i)).*(H_ssb_echoes* [beta_ue(1).*precoderSSB_GoB(:,i); beta_ue(2).*precoderSSB_GoB(:,i); beta_ue(3).*precoderSSB_GoB(:,i); beta_ue(4).*precoderSSB_GoB(:,i); beta_ue(5).*precoderSSB_GoB(:,i)] ) .*SSB_symbols1(i);
        
    end
    Y_APr(:,i) = nom + denom + noise_apr(:,i);
    y_apr_comb(i) = V(:,i)'*Y_APr(:,i);
    SINR_s(i) = abs(V(:,i)'*nom)/abs(V(:,i)'*(denom + noise_apr(:,i)) );
end

% Rx signal at UE % raw
y_ue = zeros(1,len_mat);
noise_ue = sqrt(noise_var/2).*( randn(1,len_mat) + 1i*randn(1,len_mat));
SINR_uj = zeros(1,len_mat);
for i=1:len_mat
    h_ju = 1/sqrt(2).*( randn(M,1) + 1i*randn(M,1));
    h_iu_2 = 1/sqrt(2).*( randn(M,1) + 1i*randn(M,1));
    h_iu_3 = 1/sqrt(2).*( randn(M,1) + 1i*randn(M,1));
    desire = sqrt(rho_com(i).*beta_j_ue(i)).*conj(h_ju)'*precoderSSB_GoB(:,i).*SSB_symbols1(i);
    undesire = sqrt(rho_com(i)).*(sqrt(beta_i_ue(i,1)).*conj(h_iu_2)'.*SSB_symbols2(i) + sqrt(beta_i_ue(i,2)).*conj(h_iu_3)'.*SSB_symbols3(i) )* precoderSSB_GoB(:,i) ...
        + (sqrt(beta_j_ue(i)).*conj(h_ju)'*precoder_Sens_APs(:,i,1) + sqrt(beta_i_ue(i,1)).*conj(h_iu_2)'*precoder_Sens_APs(:,i,2) + sqrt(beta_i_ue(i,2)).*conj(h_iu_3)'*precoder_Sens_APs(:,i,3) ).*Sens_symbols(i);
    y_ue(i) = desire + undesire + noise_ue(i);
    SINR_uj(i) = abs(desire)/abs(undesire + noise_ue(i));
end


% SINR with close form for opt
wj_opt = ones(M,3);
SINR_s_opt =zeros(len_mat,1); 

for i=1:1%len_mat
    nom_SINR_s = rcs_var*abs( V(:,i)'*(sqrt( beta_j_Pq(1,i)).*H1(:,:,i)*wj_opt(:,1) + sqrt(beta_j_Pq(2,i)).*H2(:,:,i)*wj_opt(:,2) + sqrt(beta_j_Pq(3,i)).*H3(:,:,i)*wj_opt(:,3) )  )^2;
    denom_SINR_s = rho_com(i).*3*sum(beta_ue) + noise_var;
    SINR_s_opt(i) = nom_SINR_s/denom_SINR_s;  % max this SENSING SINR
    
    nom_SINR_u = rho_com(i)*beta_j_ue(i);
    denom_SINR_u = rho_com(i)*sum(beta_i_ue(i,:)) + (beta_j_ue(i)*norm(wj_opt(:,1))^2 + beta_i_ue(i,1)*norm(wj_opt(:,2))^2 + beta_i_ue(i,2)*norm(wj_opt(:,3))^2) + noise_var;
    SINR_u = nom_SINR_u/denom_SINR_u;

    % such that
    tot_pwr = [norm(wj_opt(:,1))^2; norm(wj_opt(:,2))^2; norm(wj_opt(:,3))^2] + rho_com(i) ; % less than Pmax
    condi2 = abs([wj_opt(:,1)';wj_opt(:,2)';wj_opt(:,3)']*precoderSSB_GoB(:,i) ) ; % equal to zero
    condi3 =  SINR_u; %SINR_u > gamma

end











%%
figure
scatter3(AP_positions(2,:),AP_positions(1,:),AP_positions(3,:),'or')
hold on
scatter3(AP_positions(2,end),AP_positions(1,end),AP_positions(3,end),'ob')
xlabel('x'),ylabel('y'),zlabel('z')

scatter3(point_positions(2,:),point_positions(1,:),point_positions(3,:),'.k')
legend('APs (tx)','APs (rx)','Points')
xlim([-20,220]),ylim([-20,220]),zlim([0,150]) 