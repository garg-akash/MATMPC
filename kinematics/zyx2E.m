function E = zyx2E(angles)
E = [0,  sin(angles(3))/cos(angles(2)),   cos(angles(3))/cos(angles(2));
    0,       cos(angles(3)),        -sin(angles(3));
    1,  sin(angles(3))*sin(angles(2))/cos(angles(2)),   cos(angles(3))*sin(angles(2))/cos(angles(2))];
end