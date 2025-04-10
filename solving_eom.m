beta = pi/2;
m1 = 1;
m2 = 1;
l = 1;

E1 = [1,0];
E2 = [0,1];

% solving for the u-coordinates
tspan = [0 10];
y0 = [0,0,0,1,1];
[t,y] = ode45(@(t,y) GAeom(t,y,beta,m1,m2,l), tspan, y0);


figure()
hold on

animation = VideoWriter('sleigh_angled.mp4', 'MPEG-4');
animation.FrameRate = 10;
open(animation);

axis([min(y(:,1))-2*l, max(y(:,1))+2*l, min(y(:,2))-2*l, max(y(:,2))+2*l])

for i = 1:length(t)

    r1 = [y(i,1),y(i,2)];
    e1 = cos(y(i,3))*E1+sin(y(i,3))*E2;
    r2 = r1+l*e1;


    particle1 = plot(r1(1),r1(2),'Marker','o','LineWidth',2,'Color','k');
    
    particle2 = plot(r2(1),r2(2),'Marker','o','LineWidth',2,'Color','b');

    rod = plot([r1(1), r2(1)], [r1(2), r2(2)],'LineWidth',2,'Color','k');


    drawnow
    pause(0.01)
    writeVideo(animation, getframe(gcf))

    delete(particle1)
    delete(particle2)
    delete(rod)


end
close(animation)





figure()
hold on
plot(y(:,1))
title('x')

figure()
hold on
plot(y(:,2))
title('y')

figure()
hold on
plot(y(:,3))
title('theta')

figure()
hold on
plot(y(:,1),y(:,2))
title('x vx y')

% % solving for the x-coordinates
% H = [cos(beta + theta)/cos(beta), -l*sin(theta)*tan(beta), -sin(theta);
%     sin(beta + theta)/cos(beta),  l*cos(theta)*tan(beta),  cos(theta);
%     0, 1, 0];



function dydt = GAeom(t,y,beta,m1,m2,l)
    dydt = zeros(5,1);
    % y1=x, y2=y, y3=theta, y4=u1, y5=u2
    % the components of the matrix H relate the derivatives of x, y, and
    % thera to u1, u2, and u3
    dydt(1) = cos(beta + y(3))/cos(beta)*y(4)+l*sin(y(3))*y(5); 
    dydt(2) = sin(beta + y(3))/cos(beta)*y(4)-l*cos(y(3))*y(5);
    dydt(3) = y(5);
    dydt(4) = (m1*y(5)*((sin(2*beta)*y(4))/2 - l*cos(beta)^2*y(5)))/(m2 + m1*cos(beta)^2);
    dydt(5) = (2*y(5)*(m1*y(4) + m2*y(4) - (l*m1*sin(2*beta)*y(5))/2))/(l*(2*m2 + 2*m1*cos(beta)^2));
end

% function dudt = GAeom(t,u,beta)
%     dudt = zeros(2,1);
%     dudt(1) = -(m1*u2*((sin(2*beta)*u1)/2 - l*cos(beta)^2*u2))/(m2 + m1*cos(beta)^2);
%     dudt(2) = -(2*u2*(m1*u1 + m2*u1 - (l*m1*sin(2*beta)*u2)/2))/(l*(2*m2 + 2*m1*cos(beta)^2));
% end