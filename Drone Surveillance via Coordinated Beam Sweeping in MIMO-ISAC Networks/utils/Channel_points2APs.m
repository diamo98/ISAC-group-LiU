function [h_ap1_2_pi, h_ap2_2_pi, h_ap3_2_pi, h_pi_2_aprx, h_APs_2_pi, h_direct_dep, h_direct_ari] = Channel_points2APs(MV,MH, azi_deg_tx2points, ele_deg_tx2points, azi_deg_points2rx, ele_deg_points2rx, azi_direct_dep, azi_direct_ari)
    h_ap1_2_pi = genURA_Response_rad(MV,MH, deg2rad(azi_deg_tx2points(1,:)), deg2rad(ele_deg_tx2points(1,:))); % ap1
    h_ap2_2_pi = genURA_Response_rad(MV,MH, deg2rad(azi_deg_tx2points(2,:)), deg2rad(ele_deg_tx2points(2,:))); % ap2
    h_ap3_2_pi = genURA_Response_rad(MV,MH, deg2rad(azi_deg_tx2points(3,:)), deg2rad(ele_deg_tx2points(3,:))); % ap3
    h_pi_2_aprx = genURA_Response_rad(MV,MH, deg2rad(azi_deg_points2rx), deg2rad(ele_deg_points2rx)); % 
    h_direct_dep = genURA_Response_rad(MV,MH, deg2rad(azi_direct_dep), zeros(1,3) ) ;
    h_direct_ari =  genURA_Response_rad(MV,MH, deg2rad(azi_direct_ari), zeros(1,3) ) ;

    [rol, col] = size(h_ap3_2_pi);
    h_APs_2_pi = zeros(rol,col,3);
    h_APs_2_pi(:,:,1) = h_ap1_2_pi;
    h_APs_2_pi(:,:,2) = h_ap2_2_pi;
    h_APs_2_pi(:,:,3) = h_ap3_2_pi;

end