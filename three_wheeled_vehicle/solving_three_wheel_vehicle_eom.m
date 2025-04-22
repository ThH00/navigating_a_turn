% clear all
close all

global l d m1 m2 m3 beta l3 l1 delta gamma Izz
% Given values (example values, replace with your own)
l = 1.0;        % total length
m1 = 1.0;       % mass 1
m2 = 1.0;       % mass 2
m3 = 1.0;       % mass 3
beta = 2*pi/3;    % beta (though not used in equations directly)
d = 2*l^2*(1+cos(2*beta));

lx = sqrt(l^2-(d/2)^2);
l3 = lx-m3*lx/(m1+m2+m3);
l1 = sqrt((lx-l3)^2+(d/2)^2);
delta = acos((l1^2+l^2-l3^2)/(2*l1*l));
gamma = acos((l1^2+l3^2-l^2)/(2*l1*l3));
Izz = m3*l3^2+(m1+m2)*l1^2;

figure(1)
hold on
set(gcf,'color','w');   % setting background color white

% initial conditions of plotted motions
u10 = 0.5;
u20 = 1;

n = 150;
tspan = linspace(0,1,n);

[t,u] = ode45(@bheom,tspan,[u10;u20;0;0;pi/3]);

xc = u(:,3);
yc = u(:,4);
theta = u(:,5);


figure()
hold on

animation = VideoWriter('three_wheel_vehicle.mp4', 'MPEG-4');
animation.FrameRate = 10;
open(animation);

axis([min(xc)-1, max(xc)+1, min(yc)-1, max(yc)+1])
box on
axis equal

plot(xc,yc,'k','linewidth',2);

% animation
for i = 1:n

    ex_plot = quiver(xc(i),yc(i),cos(theta(i)),sin(theta(i)),'r','LineWidth',2);
    ey_plot = quiver(xc(i),yc(i),-sin(theta(i)),cos(theta(i)),'b','LineWidth',2);

    axis equal
    axis([min(xc)-1, max(xc)+1, min(yc)-1, max(yc)+1])

    drawnow
    writeVideo(animation, getframe(gcf))


    delete(ex_plot)
    delete(ey_plot)

end

close(animation)


function duxdt = bheom(t,u)
    % u1 = u(1);
    % u2 = u(2);
    % x = u(3);
    % y = u(4);
    % theta = u(5);
    global l d m1 m2 m3 beta l3 l1 delta gamma Izz
    m = m1+m2+m3;
    duxdt = [(l*cos(gamma)*u(2)*(m*cos(gamma)*u(2)*l^2 - m*sin(gamma)*u(1)*l + Izz*cos(gamma)*u(2)))/(m*cos(gamma)^2*l^2 + Izz);
        -(l*m*u(2)*(2*u(1) - l*sin(2*gamma)*u(2)))/(2*(m*cos(gamma)^2*l^2 + Izz));
        cos(gamma + u(5))/cos(gamma)*u(1)+l*sin(u(5))*u(2);
        sin(gamma + u(5))/cos(gamma)*u(1)-l*cos(u(5))*u(2);
        u(2)];
end

% H = 
% 
% [cos(gamma + theta(t))/cos(gamma),  l*sin(theta(t)), -sin(theta(t))]
% [sin(gamma + theta(t))/cos(gamma), -l*cos(theta(t)),  cos(theta(t))]
% [                               0,                1,              0]
% 
% 
% % udot = 
% (l*cos(gamma)*u(2)*(m*cos(gamma)*u(2)*l^2 - m*sin(gamma)*u(1)*l + I_zz*cos(gamma)*u(2)))/(m*cos(gamma)^2*l^2 + I_zz)
%                                            -(l*m*u(2)*(2*u(1) - l*sin(2*gamma)*u(2)))/(2*(m*cos(gamma)^2*l^2 + I_zz))