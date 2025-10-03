function [azimu_beams,polar_beams,BeamList,nr_tot_beams,precoder] = orthogonalGoB_PrecodeGen(MV)
    
    M = MV^2;
    n_arr = 1:5;
    angles = [flip(-1.*asind(2.*n_arr./MV)), 0, asind(2.*n_arr./MV)];
    
    angles_des = sort(angles,'descend');
    nr_unique_beams = length(angles);
    nr_tot_beams = 64;
    azimu_beams = angles_des;
    polar_beams = angles_des(6:end); % MODIFY ME
    BeamList = zeros(nr_unique_beams^2,2);
    
    for j=1:nr_unique_beams
        
        for i=1:nr_unique_beams
            BeamList(i+((j-1)*nr_unique_beams),:) = [angles_des(j), angles_des(i)];
            
        end
    end
    
    BeamList = BeamList(( nr_unique_beams^2-nr_tot_beams+1):end,:);
    precoder = zeros(M,nr_tot_beams);

    for i=1:nr_tot_beams
        precoder(:,i) = genURA_Response_rad(MV,MV,deg2rad(BeamList(i,2)),deg2rad(BeamList(i,1))); % azi, elevation
    end
    precoder = repelem(precoder, 1, 4);
    precoder = 1./sqrt(M).*precoder;

end