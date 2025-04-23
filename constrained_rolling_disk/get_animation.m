function [] = get_animation(q)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


load('outputs/q.mat')
ntime = size(q,1);

R = 0.3750;
h = 0;

berkeley_blue = [0, 50, 98]/256;
california_gold = [253, 181, 21]/256;
soybean = [157, 173, 51]/256;
lawrence = [0,176,218]/256;
metallic_gold = [188,155,106]/256;
bayfog = [194,185,167]/256;

x1 = q(:,1);
x2 = q(:,2);
x3 = q(:,3);
psi = q(:,4);
theta = q(:,5);
phi = q(:,6);

figure()

axis equal
view_angle = linspace(-75+180,180,ntime);
view(-45,45)
xlim([-2 2])
ylim([-2 2])
zlim([0 2])
box on

% set(gca,'XTick',[], 'YTick', [], 'ZTick', [])

animation = VideoWriter('rolling_disk.mp4', 'MPEG-4');
animation.FrameRate = 100;
open(animation);

for i = 1:4:ntime
    hold on
    view(view_angle(i),45)
    
    [disk,sides,pointA,pointP,helix,~] = ...
        plot_disk(x1(i),x2(i),x3(i),psi(i),theta(i),phi(i),berkeley_blue,california_gold);
    drawnow
    writeVideo(animation, getframe(gcf));
%     pause(0.00001)
    delete(disk)
    delete(sides)
    delete(pointA)
    delete(helix)
    
end

close(animation)

end