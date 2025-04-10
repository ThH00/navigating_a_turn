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
v2 = v1+l*thetadot*e2;

% calculating v1dot and v2dot
v1dot = diff(v1,t);
v2dot = diff(v2,t);


%%
syms u1dot u2dot real
% manually replacing the derivatives in the previous expressions by u1dot
% and u2dot

v1dot = [(cos(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) + l*cos(theta)*u2*u2, (sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + l*sin(theta)*u2*u2]

v2dot = [(cos(beta + theta)*u1dot)/cos(beta) - (sin(beta + theta)*u1*u2)/cos(beta), (sin(beta + theta)*u1dot)/cos(beta) + (cos(beta + theta)*u1*u2)/cos(beta)]
 


v1dot_squared = ((cos(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) + l*cos(theta)*u2*u2)*((cos(beta + theta)*u1dot)/cos(beta) + l*sin(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) + l*cos(theta)*u2*u2)+((sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + l*sin(theta)*u2*u2)*((sin(beta + theta)*u1dot)/cos(beta) - l*cos(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + l*sin(theta)*u2*u2);
v2dot_squared = ((cos(beta + theta)*u1dot)/cos(beta) - (sin(beta + theta)*u1*u2)/cos(beta))*((cos(beta + theta)*u1dot)/cos(beta) - (sin(beta + theta)*u1*u2)/cos(beta))+((sin(beta + theta)*u1dot)/cos(beta) + (cos(beta + theta)*u1*u2)/cos(beta))*((sin(beta + theta)*u1dot)/cos(beta) + (cos(beta + theta)*u1*u2)/cos(beta));

syms m1 m2 real
% calculate the dot product manually
S = m1/2*v1dot_squared+m2/2*v2dot_squared;

simplify(diff(S,u1dot))
simplify(diff(S,u2dot))

dSdu1dot = (m1*u1dot + m2*u1dot + l*m1*cos(beta)^2*u2^2 - (l*m1*u2dot*sin(2*beta))/2)/cos(beta)^2;
dSdu2dot = -(l*m1*(u1dot*sin(beta) - l*u2dot*cos(beta) + cos(beta)*u1*u2))/cos(beta);

A = [diff(dSdu1dot, u1dot), diff(dSdu1dot, u2dot);
    diff(dSdu2dot, u1dot), diff(dSdu2dot, u2dot)];


 
% A(t) =
% 
% [      (m1 + m2)/cos(beta)^2, -(l*m1*sin(2*beta))/(2*cos(beta)^2)]
% [-(l*m1*sin(beta))/cos(beta),                              l^2*m1]

% B = simplify([dSdu1dot;dSdu2dot]+A*[u1dot; u2dot])


% A\B
% 
% simplify(ans)
% 
