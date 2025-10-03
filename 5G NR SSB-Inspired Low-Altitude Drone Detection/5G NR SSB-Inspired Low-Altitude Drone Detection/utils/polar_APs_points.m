function [azi_deg, ele_deg, radius_n,azi_direct_dep,azi_direct_ari] = polar_APs_points(AP_positions,point_positions) % polar [azi_d; ele_d; radius]
    [~,num_tx_ap] = size(AP_positions);
    [~,num_points]=size(point_positions);
    num_tx_ap = num_tx_ap-1; % exclude rx ap

    azi_deg = zeros(num_tx_ap,num_points);
    ele_deg = zeros(num_tx_ap,num_points);
    radius_n = zeros(num_tx_ap,num_points);

    for i=1:num_tx_ap % iterate over AP-j
        azi_deg(i,:) = atan2d((point_positions(2,:) - AP_positions(2,i)), (point_positions(1,:) - AP_positions(1,i)));
        radius_n(i,:) = sqrt((point_positions(1,:) - AP_positions(1,i)).^2 + (point_positions(2,:) - AP_positions(2,i)).^2 + (point_positions(3,:) - AP_positions(3,i)).^2);
        ele_deg(i,:) = asind((point_positions(3,:) - AP_positions(3,i)) / radius_n(i,:));
    end

    azi_direct_dep = atan2d((AP_positions(2,end) - AP_positions(2,1:end-1)), (AP_positions(1,end) - AP_positions(1,1:end-1)) );
    azi_direct_ari = atan2d((AP_positions(2,1:end-1) - AP_positions(2,end)), (AP_positions(1,1:end-1) - AP_positions(1,end)) );
    

end
