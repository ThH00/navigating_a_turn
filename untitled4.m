v0 = 0.55;
g = 9.81;
R = 0.9;

vB = sqrt(v0^2+2*g*R*cosd(30));

thetadotB = sqrt(v0^2/(R*(1.15-sind(30))*0.15*R));

angle = acos(0.15*R*thetadotB/vB)

angle_deg = angle*180/pi

% 52.9


