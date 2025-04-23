% clear all
close all

E1 = [1;0];
E2 = [0;1];

global AD AB mA mB mD detla CD AC CAD ACD Izz
% Given values (example values, replace with your own)
AD = 1.0;        % total length
mA = 1.0;       % mass A
mB = 1.0;       % mass B
mD = 1.0;       % mass D
detla = pi/3;    % delta, base angle of isosceles triangular RB
AB = sqrt(2*AD^2*(1+cos(2*detla)));   % using cosine law on triangle ABD

ED = sqrt(AD^2-(AB/2)^2);
EC = mD*ED/(mA+mB+mD);
CD = ED-EC;
AC = sqrt(EC^2+(AB/2)^2);
CAD = acos((AC^2+AD^2-CD^2)/(2*AC*AD)); % angle
CBD = CAD;  % angle
ACD = acos((AC^2+CD^2-AD^2)/(2*AC*CD)); % angle
Izz = mD*CD^2+(mA+mB)*AC^2;

ACE = pi-ACD;

global beta real
beta = pi/40;

figure(1)
hold on
set(gcf,'color','w');   % setting background color white
box on

% initial conditions of plotted motions
u10 = 2;
u20 = 1;

n = 150;
tspan = linspace(0,1,n);

[t,u] = ode45(@bheom,tspan,[u10;u20;0;0;pi/3]);

xC = u(:,3);
yC = u(:,4);
theta = u(:,5);


figure()
hold on
box on

animation = VideoWriter('three_wheel_vehicle.mp4', 'MPEG-4');
animation.FrameRate = 10;
open(animation);

axis([min(xC)-1, max(xC)+1, min(yC)-1, max(yC)+1])
box on
axis equal

plot(xC,yC,'k','linewidth',2);

% animation
for i = 1:n

    e1 = cos(theta(i))*E1+sin(theta(i))*E2;
    e2 = -sin(theta(i))*E1+cos(theta(i))*E2;

    rC = xC(i)*E1+yC(i)*E2;
    rA = rC+AC*(-cos(ACE)*e1+sin(ACE)*e2);
    rB = rC+AC*(-cos(ACE)*e1-sin(ACE)*e2);
    rD = rC+CD*e1;

    AD_plot = plot([rA(1), rD(1)], [rA(2), rD(2)],'k','LineWidth',2);
    AB_plot = plot([rA(1), rB(1)], [rA(2), rB(2)],'k','LineWidth',2);
    BD_plot = plot([rD(1), rB(1)], [rD(2), rB(2)],'k','LineWidth',2);


    ex_plot = quiver(xC(i),yC(i),cos(theta(i)),sin(theta(i)),'r','LineWidth',2);
    ey_plot = quiver(xC(i),yC(i),-sin(theta(i)),cos(theta(i)),'b','LineWidth',2);

    axis equal
    axis([min(xC)-1, max(xC)+1, min(yC)-1, max(yC)+1])

    drawnow
    writeVideo(animation, getframe(gcf))


    delete(ex_plot)
    delete(ey_plot)

    delete(AD_plot)
    delete(AB_plot)
    delete(BD_plot)

end

close(animation)


function duxdt = bheom(t,u)
    % u1 = u(1);
    % u2 = u(2);
    % x = u(3);
    % y = u(4);
    % theta = u(5);
    global mA mB mD CD Izz
    global beta
    m = mA+mB+mD;
    duxdt = [(CD*cos(beta)*u(2)*(m*cos(beta)*u(2)*CD^2 - m*sin(beta)*u(1)*CD + Izz*cos(beta)*u(2)))/(m*cos(beta)^2*CD^2 + Izz);
        -(CD*m*u(2)*(2*u(1) - CD*sin(2*beta)*u(2)))/(2*(m*cos(beta)^2*CD^2 + Izz));
        cos(beta + u(5))/cos(beta)*u(1)+CD*sin(u(5))*u(2);
        sin(beta + u(5))/cos(beta)*u(1)-CD*cos(u(5))*u(2);
        u(2)];
end

% H = 
% 
% [cos(beta + theta(t))/cos(beta),  l*sin(theta(t)), -sin(theta(t))]
% [sin(beta + theta(t))/cos(beta), -l*cos(theta(t)),  cos(theta(t))]
% [                               0,                1,              0]
% 
% 
% % udot = 
% (l*cos(gamma)*u(2)*(m*cos(gamma)*u(2)*l^2 - m*sin(gamma)*u(1)*l + I_zz*cos(gamma)*u(2)))/(m*cos(gamma)^2*l^2 + I_zz)
%                                            -(l*m*u(2)*(2*u(1) - l*sin(2*gamma)*u(2)))/(2*(m*cos(gamma)^2*l^2 + I_zz))