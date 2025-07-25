function [xp,yp,zp] = mapc2m(xc,yc)

R = 1;
r = 0.4;

parms = read_vars();
lat = parms.latitude;
lng = parms.longitude;

maps = {'cubedsphere','latlong','pillowsphere'};
map = maps{parms.mapping};


switch map
    case 'cubedsphere'
        [xp,yp,zp] = mapc2m_cubedsphere(xc,yc);
    case 'latlong'
        s = 0.0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);

        % Map into [0,1]x[0,1]
        xc2 = lng(1) + (lng(2) - lng(1))*xc1;
        yc2 = lat(1) + (lat(2) - lat(1))*yc1;
        [xp,yp,zp] = mapc2m_latlong(xc2,yc2);
    case 'pillowsphere'
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);

end

% R = rotate_map(parms.theta,parms.phi);
% [m,n] = size(xp);
% 
% XYZ = R*[xp(:)'; yp(:)'; zp(:)'];
% 
% 
% 
% xp = reshape(XYZ(1,:),m,n);
% yp = reshape(XYZ(2,:),m,n);
% zp = reshape(XYZ(3,:),m,n);


end
