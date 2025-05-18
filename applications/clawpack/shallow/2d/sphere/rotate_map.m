function R = rotate_map(phi,theta)

R = eye(3);

R1 = [cos(phi), -sin(phi), 0;
    sin(phi), cos(phi), 0;
    0, 0, 1];

R2 = [cos(theta), 0, -sin(theta);
    0, 1, 0;
    sin(theta), 0, cos(theta)];

R = R2*R1;


