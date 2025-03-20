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

v1dot = diff(v1,t);
v2dot = diff(v2,t);


%%
syms u1dot u2dot real
v1dot = [(cos(beta + theta)*u1dot)/cos(beta) - (sin(beta + theta)*u1*u2)/cos(beta) - l*sin(theta)*tan(beta)*u2dot - l*cos(theta)*tan(beta)*u2*u2, (sin(beta + theta)*u1dot)/cos(beta) + (cos(beta + theta)*u1*thetadot)/cos(beta) + l*cos(theta)*tan(beta)*u2dot - l*sin(theta)*tan(beta)*u2*u2];
v2dot = [(cos(beta + theta)*u1dot)/cos(beta) + l*cos(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) - l*sin(theta)*tan(beta)*u2dot - l*sin(theta)*u2*u2 - l*cos(theta)*tan(beta)*u2*thetadot, (sin(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + l*cos(theta)*tan(beta)*u2dot + l*cos(theta)*u2*u2 - l*sin(theta)*tan(beta)*u2*u2];
 
syms m1 m2 real
% calculate the dot product manually
S = m1/2*dot(v1,v1)+m2/2*dot(v2,v2);

diff(S,u1dot)
diff(S,u2dot)