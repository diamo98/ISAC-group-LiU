function [beta_AD, beta_AU, beta_AnR, beta_ADnR ] = modify_scheme(N,frac_dis)

%M=100; % ISAC MIMO antennas
%N=50; % number of repeaters
theta_drone = pi/6; % drone angle rad
%theta_drone_deg = rad2deg(theta_drone); % drone angle deg
lambda = 3e8/15e9; rcs= 10^(-10/10); 
% function
%a = @(theta) [exp(1i*pi*(0:M-1)'*cos(theta))]; % array response
onelink_beta = @(l1) (lambda^2./(4*pi.*l1).^2);
twolink_beta = @(l1,l2) (rcs*lambda^2./((4*pi)^3.*l1.^2.*l2.^2));
% position
%frac_dis = 100;
l_AD = 500/frac_dis ; l_A1 = 250/frac_dis ; l_AU = 100/frac_dis ; 
%d = 10/frac_dis;
d = 400/N/frac_dis ; % repeater spacing
AP_pos = [0;0]; D_pos = [l_AD*cos(theta_drone) ; l_AD*sin(theta_drone)]; Repeater_pos = [l_A1+(0:N-1)*d;zeros(1,N)];
l_AnR = l_A1+(0:N-1)*d; l_DnR = vecnorm(Repeater_pos-D_pos);
% beta
beta_AD= twolink_beta(l_AD,l_AD); beta_AU= onelink_beta(l_AU); beta_AnR= onelink_beta(l_AnR)'; beta_ADnR = twolink_beta(l_AD,l_DnR)';

end



