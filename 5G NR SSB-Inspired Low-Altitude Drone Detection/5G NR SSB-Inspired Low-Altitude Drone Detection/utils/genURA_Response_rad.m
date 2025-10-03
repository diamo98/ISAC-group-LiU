function [URA_response] = genURA_Response_rad(MV,MH,azimuth,elevation)
    
    lv=1:MV;
    lh=1:MH;
    M = MV*MH;
    nr_resp = length(azimuth);
    URA_response = zeros(M,nr_resp);

    for i=1:nr_resp
        respA = [exp(-1i*pi*lv'.*sin(elevation(i)))];
        respB = [exp(-1i*pi*lh'.*sin(azimuth(i)).*cos(elevation(i)))];
        URA_response(:,i) = kron(respA,respB);
    end

end