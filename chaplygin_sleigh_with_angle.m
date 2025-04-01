% obtaining the Gibbs-Appell equations of motion

% finding the matrix H such that
% [xdot; ydot; thetadot] = H*[u1; u2; u3]

syms l beta real
syms theta(t)
syms xdot ydot thetadot real

E1 = [1,0];
E2 = [0,1];

e1 = cos(theta)*E1+sin(theta)*E2;
e2 = -sin(theta)*E1+cos(theta)*E2;

v1 = xdot*E1+ydot*E2;
v2 = v1+l*thetadot*e2;

u1 = xdot*cos(theta)+ydot*sin(theta); %dot(v1,e1);
u2 = thetadot;
u3 = (xdot - l*thetadot*sin(theta))*(- sin(theta) - cos(theta)*tan(beta))+(ydot + l*thetadot*cos(theta))*(cos(theta) - sin(theta)*tan(beta)); %dot(v2,e2-e1*tan(beta));

J = [diff(u1,xdot), diff(u1, ydot), diff(u1,thetadot);
    diff(u2,xdot), diff(u2, ydot), diff(u2,thetadot);
    diff(u3,xdot), diff(u3, ydot), diff(u3,thetadot)];

H = simplify(inv(J));


% writing v1 and v2 in terms of u
syms t
syms u1(t)
syms u2(t)
syms u3(t)

xdot = cos(beta + theta(t))/cos(beta)*u1+l*sin(theta(t))*u2; %H(1,:)*[u1;u2;0];
ydot = sin(beta + theta(t))/cos(beta)*u1-l*cos(theta(t))*u2; %H(2,:)*[u1;u2;0];
thetadot = u2; %H(3,:)*[u1;u2;0];

v1 = xdot*E1+ydot*E2;
v2 = v1+l*thetadot*e1;

% calculating v1dot and v2dot
v1dot = diff(v1,t);
v2dot = diff(v2,t);


%%
syms u1dot u2dot real
% manually replacing the derivatives in the previous expressions by u1dot
% and u2dot
 
v1dot = [(cos(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot - (sin(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot, (sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*sin(theta)*u2*thetadot];
 
v2dot = [(cos(beta + theta)*u1dot)/cos(beta) + l*cos(theta)*u2dot + l*sin(theta)*u2dot - (sin(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot - l*sin(theta)*u2*thetadot, (sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + l*sin(theta)*u2dot + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot + l*sin(theta)*u2*thetadot];

v1dot_squared = ((cos(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot - (sin(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot)*((cos(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot - (sin(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot)+((sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*sin(theta)*u2*thetadot)*((sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*sin(theta)*u2*thetadot);
v2dot_squared = ((cos(beta + theta)*u1dot)/cos(beta) + l*cos(theta)*u2dot + l*sin(theta)*u2dot - (sin(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot - l*sin(theta)*u2*thetadot)*((cos(beta + theta)*u1dot)/cos(beta) + l*cos(theta)*u2dot + l*sin(theta)*u2dot - (sin(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot - l*sin(theta)*u2*thetadot)+((sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + l*sin(theta)*u2dot + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot + l*sin(theta)*u2*thetadot)*((sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + l*sin(theta)*u2dot + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*u2*thetadot + l*sin(theta)*u2*thetadot);



syms m1 m2 real
% calculate the dot product manually
S = m1/2*v1dot_squared+m2/2*v2dot_squared;

simplify(diff(S,u1dot))
simplify(diff(S,u2dot))



dSdu1dot = (m1*u1dot + m2*u1dot + l*m1*cos(beta)^2*u2^2 + l*m2*cos(beta)^2*u2^2 + (l*m2*sin(2*beta)*u2^2)/2 + l*m2*u2dot*cos(beta)^2 - (l*m1*u2dot*sin(2*beta))/2 - (l*m2*u2dot*sin(2*beta))/2)/cos(beta)^2;
 
dSdu2dot = -(l*m1*u1dot*sin(beta) - l*m2*u1dot*cos(beta) + l*m2*u1dot*sin(beta) - l^2*m1*u2dot*cos(beta) - 2*l^2*m2*u2dot*cos(beta) + l*m1*cos(beta)*u1(t)*u2(t) + l*m2*cos(beta)*u1(t)*u2(t) + l*m2*sin(beta)*u1(t)*u2(t))/cos(beta);


A = [diff(dSdu1dot, u1dot), diff(dSdu1dot, u2dot);
    diff(dSdu2dot, u1dot), diff(dSdu2dot, u2dot)];

 
% simplify(A)
% 
% ans(t) =
% 
% [                                      (m1 + m2)/cos(beta)^2, -(l*(m1*sin(beta) - m2*cos(beta) + m2*sin(beta)))/cos(beta)]
% [-(l*(m1*sin(beta) - m2*cos(beta) + m2*sin(beta)))/cos(beta),                                             l^2*(m1 + 2*m2)]
% 


% B = simplify([dSdu1dot;dSdu2dot]-A*[u1dot; u2dot])
% ans(t) =
% 
%      (l*u2(t)^2*(m1*cos(beta) + m2*cos(beta) + m2*sin(beta)))/cos(beta)
% -(l*u1(t)*u2(t)*(m1*cos(beta) + m2*cos(beta) + m2*sin(beta)))/cos(beta)

B = [(l*u2(t)^2*(m1*cos(beta) + m2*cos(beta) + m2*sin(beta)))/cos(beta); -(l*u1(t)*u2(t)*(m1*cos(beta) + m2*cos(beta) + m2*sin(beta)))/cos(beta)];
% 
% A\B
% 
% ans(t) =
% 
%              (u2(t)*(m1*cos(beta) + m2*cos(beta) + m2*sin(beta))*(2*m2*cos(beta)^2*u1(t) - m1*sin(2*beta)*u1(t) - m2*sin(2*beta)*u1(t) + 2*l*m1*cos(beta)^2*u2(t) + 4*l*m2*cos(beta)^2*u2(t)))/(2*m1^2*cos(beta) - 2*m2^2*cos(beta)^3 + 4*m2^2*cos(beta) + m2^2*sin(2*beta)*cos(beta) + 2*m2^2*cos(beta)^2*sin(beta) - m1^2*sin(2*beta)*sin(beta) - m2^2*sin(2*beta)*sin(beta) + 6*m1*m2*cos(beta) + m1*m2*sin(2*beta)*cos(beta) + 2*m1*m2*cos(beta)^2*sin(beta) - 2*m1*m2*sin(2*beta)*sin(beta))
% -(2*u2(t)*(m1*cos(beta) + m2*cos(beta) + m2*sin(beta))*(m1*u1(t) + m2*u1(t) + l*m2*cos(beta)^2*u2(t) - l*m1*cos(beta)*sin(beta)*u2(t) - l*m2*cos(beta)*sin(beta)*u2(t)))/(2*l*m1^2*cos(beta) + 4*l*m2^2*cos(beta) - 2*l*m2^2*cos(beta)^3 + 6*l*m1*m2*cos(beta) + l*m2^2*sin(2*beta)*cos(beta) + 2*l*m2^2*cos(beta)^2*sin(beta) - l*m1^2*sin(2*beta)*sin(beta) - l*m2^2*sin(2*beta)*sin(beta) + l*m1*m2*sin(2*beta)*cos(beta) + 2*l*m1*m2*cos(beta)^2*sin(beta) - 2*l*m1*m2*sin(2*beta)*sin(beta))
% 
% 
% function dudt = GAeom(t,u,beta)
%     dudt = zeros(2,1);
%     dudt(1) = (u2*(m1*cos(beta) + m2*cos(beta) + m2*sin(beta))*(2*m2*cos(beta)^2*u1 - m1*sin(2*beta)*u1 - m2*sin(2*beta)*u1 + 2*l*m1*cos(beta)^2*u2 + 4*l*m2*cos(beta)^2*u2))/(2*m1^2*cos(beta) - 2*m2^2*cos(beta)^3 + 4*m2^2*cos(beta) + m2^2*sin(2*beta)*cos(beta) + 2*m2^2*cos(beta)^2*sin(beta) - m1^2*sin(2*beta)*sin(beta) - m2^2*sin(2*beta)*sin(beta) + 6*m1*m2*cos(beta) + m1*m2*sin(2*beta)*cos(beta) + 2*m1*m2*cos(beta)^2*sin(beta) - 2*m1*m2*sin(2*beta)*sin(beta));
%     dudt(2) = -(2*u2*(m1*cos(beta) + m2*cos(beta) + m2*sin(beta))*(m1*u1 + m2*u1 + l*m2*cos(beta)^2*u2 - l*m1*cos(beta)*sin(beta)*u2 - l*m2*cos(beta)*sin(beta)*u2))/(2*l*m1^2*cos(beta) + 4*l*m2^2*cos(beta) - 2*l*m2^2*cos(beta)^3 + 6*l*m1*m2*cos(beta) + l*m2^2*sin(2*beta)*cos(beta) + 2*l*m2^2*cos(beta)^2*sin(beta) - l*m1^2*sin(2*beta)*sin(beta) - l*m2^2*sin(2*beta)*sin(beta) + l*m1*m2*sin(2*beta)*cos(beta) + 2*l*m1*m2*cos(beta)^2*sin(beta) - 2*l*m1*m2*sin(2*beta)*sin(beta));
% end