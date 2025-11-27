function [X,Y,Z] = polar2cartisian(azimuth,polar,radius)
    X = radius.*sin(polar).*cos(azimuth);
    Y = radius.*sin(polar).*sin(azimuth);
    Z = radius.*(cos(polar));
end