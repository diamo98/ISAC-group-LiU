function [AP_positions,point_positions] = setup_cart(r,z,point_row,point_col,point_deep,d,plot_flag)
    %R = r/cosd(30);
    AP1 = [0,0,0]';
    AP2 = [0,2*r,0]';
    APr = [3*r/(2*cosd(30)),r,0]'; 
    AP3 = [3*r/(2*cosd(30)),3*r,0]';
    AP_positions = [AP1,AP2,AP3,APr];
    point_positions = zeros(3,point_row*point_col*point_deep);
    count = 1;
    x=r; y=3/2*r;
    for i=1:point_deep
        for j=1:point_row
            for k=1:point_col
                point_positions (:,count) = [x-(d*(k-1)); y-(d*(j-1)); z+(d*(i-1))];
                count = count + 1;
            end
        end
    end
    %plot_flag=1;
    if plot_flag==1
        figure
        scatter3(AP_positions(1,:),AP_positions(2,:),AP_positions(3,:),'or')
        hold on
        scatter3(AP_positions(1,end),AP_positions(2,end),AP_positions(3,end),'ob')
        xlabel('x'),ylabel('y'),zlabel('z')
        
        scatter3(point_positions(1,:),point_positions(2,:),point_positions(3,:),'.k')
        legend('AP-j','AP-r','Points')
        xlim([-20,max(AP_positions(1,:)+20)]),ylim([-20,max(AP_positions(2,:))+20]),zlim([0,50]) 
    end
end