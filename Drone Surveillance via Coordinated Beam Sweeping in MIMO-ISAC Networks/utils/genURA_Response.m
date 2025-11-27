function [URA_response] = genURA_Response(MV,MH,azimuth,elevation)
    
    lv=1:MV;
    lh=1:MH;
    respA = [exp(-1i*pi*lv'.*sin(elevation))];
    respB = [exp(-1i*pi*lh'.*sin(azimuth).*cos(elevation))];
    URA_response = kron(respA,respB);

    

end