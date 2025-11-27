%% 2D gain
MV=10;
MH=10;
[azimu_beams,polar_beams,BeamList,nr_tot_beams,precoder] = orthogonalGoB_PrecodeGen(MV);
%
%[X,Y,Z,Gain_mat_abs] = genBeamSSB(BeamList(1,1),BeamList(1,2));
azimuth = deg2rad(BeamList(:,2));
elevation = deg2rad(BeamList(:,1));
[URA_response] = genURA_Response_multi(10,10,azimuth,elevation);

resol_1=1*1e3;
resol_2=1*1e3;
angles_azi = linspace(-pi/3,pi/3,resol_1); % all possible angles
angles_ele = deg2rad(60)*ones(1,resol_2);%pi/3*ones(1,resol_2);%linspace(-pi/2,pi/2,resol_2); % at what angle we consider
% 53.1301
angles_azi_deg = rad2deg(angles_azi);

num=9;
gain_array_1 = zeros(resol_2,num);
arbitary_resp = genURA_Response_multi(10,10,angles_azi,angles_ele);

for j=1:9
    
    
    gain_array_1(:,j) = abs(arbitary_resp'*URA_response(:,j));
    
end

figure('Position', [100, 100, 360, 300])
Markers = {'k-','k:','k--','b-','b:','b--','r-','r:','r--'};
for i=1:9
%plot(angles_azi_deg,gain_array_1),legend('Beam 1','Beam 2','Beam 3','Beam 4','Beam 5','Beam 6','Beam 7','Beam 8','Beam 9'),ylabel('Beamforming gain'),xlabel('Sweep degree in azimuth')
plot(angles_azi_deg,10*log10(gain_array_1(:,i)),strcat(Markers{i}),LineWidth=1.2)%,legend('Beam 0','Beam 2','Beam 3','Beam 4','Beam 5','Beam 6','Beam 7','Beam 8','Beam 9')
ylabel('Beamforming gain [dB]','FontSize',14),xlabel('Azimuth [degrees]','FontSize',14)
%title('10x10 UPA gain pattern at elevation = 53 degree')
hold on
end
ylim([10,20])
grid on

%% 2D gain


[azimu_beams,polar_beams,BeamList,nr_tot_beams,precoder] = orthogonalGoB_PrecodeGen(MV,MH);

%[X,Y,Z,Gain_mat_abs] = genBeamSSB(BeamList(1,1),BeamList(1,2));
azimuth = deg2rad(BeamList(:,2));
elevation = deg2rad(BeamList(:,1));
[URA_response] = genURA_Response_multi(10,10,azimuth,elevation);

resol_1=1*1e3;
resol_2=1*1e3;
angles_azi = linspace(-pi/3,pi/3,resol_1); % all possible angles
angles_ele = 0*ones(1,resol_2);%linspace(-pi/2,pi/2,resol_2); % at what angle we consider
angles_azi_deg = rad2deg(angles_azi);

num=9;
gain_array_1 = zeros(resol_2,num);
arbitary_resp = genURA_Response_multi(10,10,angles_azi,angles_ele);

for j=1:9
    
    
    gain_array_1(:,j) = abs(arbitary_resp'*URA_response(:,j+36));
    
end

figure('Position', [100, 100, 360, 300])
Markers = {'k-','k:','k--','b-','b:','b--','r-','r:','r--'};
for i=1:9
%plot(angles_azi_deg,gain_array_1),legend('Beam 1','Beam 2','Beam 3','Beam 4','Beam 5','Beam 6','Beam 7','Beam 8','Beam 9'),ylabel('Beamforming gain'),xlabel('Sweep degree in azimuth')
plot(angles_azi_deg,10*log10(gain_array_1(:,i)),strcat(Markers{i}),LineWidth=1.2)%,legend('Beam 0','Beam 2','Beam 3','Beam 4','Beam 5','Beam 6','Beam 7','Beam 8','Beam 9'),ylabel('Beamforming gain [dB]'),xlabel('Sweep degree in azimuth')
%title('10x10 UPA gain pattern at elevation = 53 degree')
hold on
end
ylim([10,20])
ylabel('Beamforming gain [dB]','FontSize',14),xlabel('Azimuth [degrees]','FontSize',14)
grid on

%%
plotflag=1;

if plotflag==1
    figure()
end

resol_1=3*1e2-1;
acc_gain = zeros(resol_1,resol_1,81);
edge_polar = pi/3;

azimuth = (BeamList(:,2));
elevation = (BeamList(:,1));

for i=1:2
    hold on
    i_in=i;%+36;
    [X,Y,Z,Gain_mat_abs] = genBeamSSB((azimuth(i_in)),(elevation(i_in)),resol_1,edge_polar,plotflag);
    colorData = X;
    surf(X,Y,Z)
    acc_gain(:,:,i_in) = Gain_mat_abs;
    i
end

%%
Min_gain = max(acc_gain,[],3);
figure
surf(Min_gain,'edgecolor','none')
%% 2D gain plot 
M = MV*MH;
grid_profile = zeros(resol_1,resol_1);
circle_gain_db = 3;
circle_gain_linear = 10^(circle_gain_db/10);
line_gain_min = M/circle_gain_linear;
azi_ele=((1:resol_1)/resol_1-0.5)*rad2deg(edge_polar)*2;
line_thresh = 1.8;

for i=1:81
    grid_profile = grid_profile+ (abs(acc_gain(:,:,i)-line_gain_min)<line_thresh);
end

figure
imagesc(azi_ele,azi_ele,grid_profile)
hold on
scatter(rad2deg(azimuth),rad2deg(elevation),'.r')
title('Grid of beams gain M - 3 dB'),xlabel('Azimuth [degrees]'),ylabel('Elevation [degrees]')


grid on


%% Gen Beamforming gain 3D UPA
function [X,Y,Z,Gain_mat_abs] = genBeamSSB(azimu_beam_in,eleva_beam_in,resol_1,edge_polar,plotflag)
    LV=10;
    LH=10;
    L=LV*LH;
    resol_1=resol_1;%4*1e2-1;
    resol_2=resol_1;%4*1e2-1;
    azimu_beam = deg2rad(azimu_beam_in);
    eleva_beam = deg2rad(eleva_beam_in);
    %m=0:(L-1);
    
    A_in = genURA_Response_rad(LV,LH,azimu_beam,eleva_beam);
    
    angles_azi = linspace(-edge_polar,edge_polar,resol_1); % all possible angles
    angles_ele = linspace(-edge_polar,edge_polar,resol_2); % at what angle we consider
    
    %gain sweep
    Gain_mat = zeros(resol_1,resol_2);
    for i=1:resol_1
        for j=1:resol_2
            Gain_mat(i,j) = A_in'*genURA_Response_rad(LV,LH,angles_azi(i),angles_ele(j));
        end
    end
    
    Gain_mat_abs = abs(Gain_mat); 
    
    [angles_ele_mesh,angles_azi_mesh] = meshgrid(angles_ele,angles_azi);
    
    angles_polar_mesh = elevation2polar(angles_ele_mesh);
    [X,Y,Z] = polar2cartisian(angles_azi_mesh,angles_polar_mesh,Gain_mat_abs);
    

    if plotflag==1
        C=(sqrt(X.^2+Y.^2+Z.^2));
        surf(X,Y,Z,C)
    
        cb = colorbar;
        cb.Label.String = 'Beamforming Gain';
        hold on
        xlabel('x')
        ylabel('y')
        zlabel('Z')
    end

end


%%


MV=10;
lv=1:MV;
respA = [exp(-1i*pi*lv'.*sin([pi/8,pi/7]))];
y = linspace(-pi/3,pi/3,1e3);
respB = [exp(-1i*pi*lv'.*sin(y))];

sum(abs(respB'*respA),2)
