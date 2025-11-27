function [beta_j_ue] = beta_j_ue_all(r,AP_positions,lambda,case_interf)
r_outer = r/cosd(30);
beta_j_ue = zeros(3,3); % row=AP, col=links
% AP-1
beta_j_ue(1,1) =  1/(r_outer)^2; % lambda^2/(4*pi*(r_outer))^2;
%r12=norm([r;0;0]-AP_positions(:,2));
if case_interf ==1 % freq reuse=3
    beta_j_ue(1,3)= 1/(2*r_outer)^2; %lambda^2/(4*pi*(2*r_outer))^2;
    beta_j_ue(1,2)= beta_j_ue(1,3); %lambda^2/(4*pi*(2*r_outer))^2;

    beta_j_ue(2,3)= beta_j_ue(1,3);
    beta_j_ue(2,1)= beta_j_ue(1,2);

    beta_j_ue(3,1)=beta_j_ue(1,3);
    beta_j_ue(3,2)=beta_j_ue(1,2);

elseif case_interf==0
    beta_j_ue(1,3)= 1/(1*r_outer)^2; %lambda^2/(4*pi*(2*r_outer))^2;
    beta_j_ue(1,2)= beta_j_ue(1,3); %lambda^2/(4*pi*(2*r_outer))^2;

    beta_j_ue(2,3)= beta_j_ue(1,3);
    beta_j_ue(2,1)= beta_j_ue(1,3);

    beta_j_ue(3,1)= beta_j_ue(1,3);
    beta_j_ue(3,2)= beta_j_ue(1,3);

end
% AP-2,3
beta_j_ue(2,2) = beta_j_ue(1,1);

beta_j_ue(3,3) = beta_j_ue(1,1);

end

% % AP-1
% beta_j_ue(1,1) = lambda^2/(4*pi*(r))^2;
% r12=norm([r;0;0]-AP_positions(:,2));
% beta_j_ue(1,3)= lambda^2/(4*pi*(r12))^2;
% r13= norm([r;0;0]-AP_positions(:,3));
% beta_j_ue(1,2)= lambda^2/(4*pi*(r13))^2;
% % AP-2,3
% beta_j_ue(2,2) = lambda^2/(4*pi*(r))^2;
% r21=norm([ISD-r;0;0]-AP_positions(:,2));
% beta_j_ue(2,3)= lambda^2/(4*pi*(r21))^2;
% r31= norm([r;0;0]-AP_positions(:,3));
% beta_j_ue(2,1)= lambda^2/(4*pi*(r31))^2;
% 
% beta_j_ue(3,3) = lambda^2/(4*pi*(r))^2;
% r21=norm([ISD-r;0;0]-AP_positions(:,2));
% beta_j_ue(3,2)= lambda^2/(4*pi*(r21))^2;
% r31= norm([r;0;0]-AP_positions(:,3));
% beta_j_ue(3,1)= lambda^2/(4*pi*(r31))^2;
% 