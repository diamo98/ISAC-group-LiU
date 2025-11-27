function [azi_deg_points2rx, ele_deg_points2rx, radius_n_points2rx] = polar_Points_APrx(AP_positions,point_positions)
    [~,num_ap] = size(AP_positions);
    [~,num_points]=size(point_positions);
    %num_tx_ap = num_tx_ap-1 % exclude rx ap
    rx_ap_index = num_ap;
    %azi_deg_points2rx = atan2d(( AP_positions(2,rx_ap_index)- point_positions(2,:) ), ( AP_positions(1,rx_ap_index) - point_positions(1,:) ));
    azi_deg_points2rx = atan2d(( point_positions(2,:)- AP_positions(2,rx_ap_index) ), ( point_positions(1,:) - AP_positions(1,rx_ap_index) ));
    % azi_deg_points2rx original one
    radius_n_points2rx = sqrt((point_positions(1,:) - AP_positions(1,rx_ap_index)).^2 + (point_positions(2,:) - AP_positions(2,rx_ap_index)).^2 + (point_positions(3,:) - AP_positions(3,rx_ap_index)).^2);
    ele_deg_points2rx = asind((point_positions(3,:) - AP_positions(3,rx_ap_index)) ./ radius_n_points2rx); % positive elevation
    %ele_deg_points2rx = asind(( AP_positions(3,rx_ap_index) -
    %point_positions(3,:) ) ./ radius_n_points2rx); % % negative elevation

end


