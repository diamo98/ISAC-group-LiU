%% Gen Beamforming gain 3D UPA
function [X,Y,Z,Gain_mat_abs] = plotBeampattern(A_in,edge_polar,plotflag)
    LV=12;
    LH=12;
    L=LV*LH;
    resol_1=4*1e2-1;
    resol_2=4*1e2-1;
    %m=0:(L-1);
    
    %A_in = genURA_Response(LV,LH,azimu_beam,eleva_beam);
    
    angles_azi = linspace(-edge_polar,edge_polar,resol_1); % all possible angles
    angles_ele = linspace(-edge_polar,edge_polar,resol_2); % at what angle we consider
    
    %gain sweep
    Gain_mat = zeros(resol_1,resol_2);
    for i=1:resol_1
        for j=1:resol_2
            Gain_mat(i,j) = A_in'*genURA_Response(LV,LH,angles_azi(i),angles_ele(j));
        end
    end
    
    Gain_mat_abs = abs(Gain_mat); 
    
    [angles_ele_mesh,angles_azi_mesh] = meshgrid(angles_ele,angles_azi);
    
    angles_polar_mesh = elevation2polar(angles_ele_mesh);
    [X,Y,Z] = polar2cartisian(angles_azi_mesh,angles_polar_mesh,Gain_mat_abs);
    

    if plotflag==1
        figure()
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
