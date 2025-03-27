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
v2 = v1+l*thetadot*e1;

u1 = xdot*cos(theta)+ydot*sin(theta); %dot(v1,e1);
u2 = thetadot;
u3 = (ydot + l*thetadot*sin(theta))*(cos(theta) - tan(beta)*sin(theta)) - (sin(theta) + tan(beta)*cos(theta))*(xdot + l*thetadot*cos(theta)); %dot(v2,e2-e1*tan(beta));

J = [diff(u1,xdot), diff(u1, ydot), diff(u1,thetadot);
    diff(u2,xdot), diff(u2, ydot), diff(u2,thetadot);
    diff(u3,xdot), diff(u3, ydot), diff(u3,thetadot)];

H = simplify(inv(J));


% writing v1 and v2 in terms of u
syms t
syms u1(t)
syms u2(t)
syms u3(t)

xdot = cos(beta + theta)/cos(beta)*u1-l*sin(theta(t))*tan(beta)*u2; %H(1,:)*[u1;u2;0];
ydot = sin(beta + theta)/cos(beta)*u1+l*cos(theta(t))*tan(beta)*u2;
thetadot = u2; %H(3,:)*[u1;u2;0];

v1 = xdot*E1+ydot*E2;
v2 = v1+l*thetadot*e1;

% calculating v1dot and v2dot
v1dot = diff(v1,t);
v2dot = diff(v2,t);


%%
syms u1dot u2dot real
v1dot = [(cos(beta + theta)*u1dot)/cos(beta) - (sin(beta + theta)*u1*u2)/cos(beta) - l*sin(theta)*tan(beta)*u2dot - l*cos(theta)*tan(beta)*u2*u2, (sin(beta + theta)*u1dot)/cos(beta) + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*tan(beta)*u2dot - l*sin(theta)*tan(beta)*u2*u2];
v2dot = [(cos(beta + theta)*u1dot)/cos(beta) + l*cos(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) - l*sin(theta)*tan(beta)*u2dot - l*sin(theta)*u2*u2 - l*cos(theta)*tan(beta)*u2*thetadot, (sin(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + l*cos(theta)*tan(beta)*u2dot + l*cos(theta)*u2*u2 - l*sin(theta)*tan(beta)*u2*u2];

v1dot_squared = ((cos(beta + theta)*u1dot)/cos(beta) - (sin(beta + theta)*u1*u2)/cos(beta) - l*sin(theta)*tan(beta)*u2dot - l*cos(theta)*tan(beta)*u2*u2)^2+((sin(beta + theta)*u1dot)/cos(beta) + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*tan(beta)*u2dot - l*sin(theta)*tan(beta)*u2*u2)^2;
v2dot_squared = ((cos(beta + theta)*u1dot)/cos(beta) + l*cos(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) - l*sin(theta)*tan(beta)*u2dot - l*sin(theta)*u2*u2 - l*cos(theta)*tan(beta)*u2*thetadot)^2+((sin(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + l*cos(theta)*tan(beta)*u2dot + l*cos(theta)*u2*u2 - l*sin(theta)*tan(beta)*u2*u2)^2;

syms m1 m2 real
% calculate the dot product manually
S = m1/2*v1dot_squared+m2/2*v2dot_squared;

simplify(diff(S,u1dot))
simplify(diff(S,u2dot))

 


eq1 = (- l*m1*cos(beta)*sin(beta)*u2(t)^2 + m1*u1dot + m2*u1dot + l*m1*u2dot + l*m2*u2dot - l*m1*u2dot*cos(beta)^2)/cos(beta)^2;
 
eq2 = (l*(2*m1*u1dot + 2*m2*u1dot - 2*m1*u1dot*cos(beta)^2 + 2*l*m1*u2dot + 2*l*m2*u2dot - 2*l*m1*u2dot*cos(beta)^2 + 2*m1*cos(beta)*sin(beta)*u1(t)*u2(t)))/(2*cos(beta)^2);

A = [diff(eq1, u1dot), diff(eq1, u2dot);
    diff(eq2, u1dot), diff(eq2, u2dot)];
% 
% A =
% 
% [                               (m1 + m2)/cos(beta)^2,               (l*m1 + l*m2 - l*m1*cos(beta)^2)/cos(beta)^2]
% [(l*(2*m1 + 2*m2 - 2*m1*cos(beta)^2))/(2*cos(beta)^2), (l*(2*l*m1 + 2*l*m2 - 2*l*m1*cos(beta)^2))/(2*cos(beta)^2)]
% 
% [eq1; eq2]-A*[u1dot;u2dot]
% 
% B=    [-(l*m1*sin(beta)*u2(t)^2)/cos(beta); (l*m1*sin(beta)*u1(t)*u2(t))/cos(beta)];
% 
% 
% A\B
% 
% ans =
% 
%                                                                                          -(sin(beta)*u2(t)*(u1(t) + l*u2(t)))/cos(beta)
% (sin(beta)*u2(t)*(m1*u1(t) + m2*u1(t) + l*m1*u2(t) + l*m2*u2(t) - l*m1*cos(beta)^2*u2(t)))/(cos(beta)*(l*m1 + l*m2 - l*m1*cos(beta)^2))
% 
% 
% function dudt = GAeom(t,u,beta)
%     dudt = zeros(2,1);
%     dudt(1) = -(sin(beta)*u2(t)*(u1(t) + l*u2(t)))/cos(beta);
%     dudt(2) = (sin(beta)*u2(t)*(m1*u1(t) + m2*u1(t) + l*m1*u2(t) + l*m2*u2(t) - l*m1*cos(beta)^2*u2(t)))/(cos(beta)*(l*m1 + l*m2 - l*m1*cos(beta)^2));
% end