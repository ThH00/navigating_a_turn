% PhaseDiagram_Movie.m
% Movie: Phase portrait and motion animation obtained from the
% Boltzmann-Hammel equations of motion
% last modified: 04/01/21
% Theresa Honein

clear all
close all

global m1 m2 ell
m1 = 0.5;
m2 = 0.5;
ell = 1;

figure(1)
hold on
set(gcf,'color','w');   % setting background color white

% plotting the phase portrait
subplot(2,2,[1,3])
hold on
axis equal
box off
axis off
axis([-0.7, 0.9, -1, 1.4]);
for i = [0,.05, .1,.15,.3,.4,.5,.6]
    for j = [-0.01,.01]
        [t,u] = ode45(@bheom,[0:0.001:100],[i;j;0;0;pi/3]);
        plot(u(:,1),ell*u(:,2),'color','k','linewidth',2)
    end
end
% axes and labels
quiver(-.7,0,1.5,0,'k','linewidth',2,'AutoScale','off')     % horizontal axis
quiver(0,-1,0,2.1,'k','linewidth',2,'AutoScale','off')      % vertical axis
text(0.7,0.1,'$u_1$','Interpreter','latex','linewidth',2,'FontSize', 24); % horizontal axis label
text(0.1,1,'$\ell u_2 = \ell \dot{\theta}$','Interpreter','latex','linewidth',2,'FontSize', 24);    % vetical axis label

% initial conditions of plotted motions
u10 = [.5,0.5];
u20 = [0.1;-0.1];

% animation
myVideo = VideoWriter('movie7');
myVideo.FrameRate = 20;
open(myVideo)
for j = 1:length(u10)
    [t,u] = ode45(@bheom,[0:0.1:12],[u10(j);u20(j);0;0;pi/3]);
    u1 = u(:,1);
    u2 = u(:,2);
    x1 = u(:,3);
    y1 = u(:,4);
    theta = u(:,5);
    x2 = x1+ell*cos(theta);
    y2 = y1+ell*sin(theta);
    
    hold on
    subplot(2,2,j*2)
    hold on
    axis('equal')
    axis off
    box off
    axis([min([x1;x2])-0.1 max([x1;x2])+0.1 min([y1;y2])-0.1 max([y1;y2])+0.1])
    % plotting trajectories
    tbh1 = plot(x1,y1,'color',[0,0+0.5*(j==1),0+1*(j==2)],'linewidth',2);
    tbh2 = plot(x2,y2,'--','color',[0,0+0.5*(j==1),0+1*(j==2)],'linewidth',2);
    % axis and labels
    quiver(-0.1,-0.1,1,0,'k','linewidth',2,'AutoScale','off')     % horizontal axis
    quiver(-0.1,-0.1,0,1,'k','linewidth',2,'AutoScale','off')      % vertical axis
    text(1.1,-0.1,'$x$','Interpreter','latex','FontSize', 24);  % horizontal axis label
    text(-0.1,1.2,'$y$','Interpreter','latex','FontSize', 24);    % vetical axis label
    for i = 2:length(t)
        subplot(2,2,j*2)
        hold on
        p1 = plot(x1(i),y1(i),'o','color',[0,0+0.5*(j==1),0+1*(j==2)],'markersize',10,...
            'MarkerFaceColor',[0,0+0.5*(j==1),0+1*(j==2)]);
        p2 = plot(x2(i),y2(i),'o','color',[0,0+0.5*(j==1),0+1*(j==2)],'markersize',10,'linewidth',2);
        sled = plot([x1(i), x2(i)],[y1(i), y2(i)],'color','k','linewidth',2);
        % leg1 = legend([tbh1 tbh2],{'trajectory of $m_1$','trajectory of $m_2$'});% ,'Location','southeast');
        % set(leg1,'Interpreter','latex');
        subplot(2,2,[1,3])
        hold on
        axis off
        box off
        c1 = plot(u1(i),u2(i),'.','markersize',5,'color',[0,0+0.5*(j==1),0+1*(j==2)]);
        c2 = plot(u1(i),u2(i),'o','markersize',10,'MarkerFaceColor',[0,0+0.5*(j==1),0+1*(j==2)]);
        pause(t(i)-t(i-1))
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
        delete(sled)
        delete(p1)
        delete(p2)
        delete(c2)
    end
end
close(myVideo)

function duxdt = bheom(t,ux)
    % u1 = ux(1);
    % u2 = ux(2);
    % x = ux(3);
    % y = ux(4);
    % theta = ux(5);
    global m1 m2 ell
    duxdt = [-m1/(m1+m2)*ell*ux(2)^2;
        ux(1)*ux(2)/ell;
        [cos(ux(5)) ell*sin(ux(5)) -sin(ux(5));
        sin(ux(5)) -ell*cos(ux(5)) cos(ux(5));
        0 1 0]*[ux(1);
        ux(2);
        0]];
end