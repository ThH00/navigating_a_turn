% obtaining the Gibbs-Appell equations of motion

%% Defining the quasi-coordinates
% finding the matrix H such that
% [xdot; ydot; thetadot] = H*[u1; u2; u3]

syms CD beta real
syms theta(t)
syms xdot ydot thetadot real

E1 = [1,0];
E2 = [0,1];

e1 = cos(theta)*E1+sin(theta)*E2;
e2 = -sin(theta)*E1+cos(theta)*E2;

vC = xdot*E1+ydot*E2;
vD = vC+CD*thetadot*e2;

u1 = xdot*cos(theta)+ydot*sin(theta); %dot(vC,e1);
u2 = thetadot;
u3 = CD*thetadot - tan(beta)*(xdot*cos(theta) + ydot*sin(theta)) + ydot*cos(theta) - xdot*sin(theta); %dot(vC,e2)+l*thetadot-tan(gamma)*dot(vC,e1);

J = [diff(u1,xdot), diff(u1, ydot), diff(u1,thetadot);
    diff(u2,xdot), diff(u2, ydot), diff(u2,thetadot);
    diff(u3,xdot), diff(u3, ydot), diff(u3,thetadot)];

H = simplify(inv(J));


% writing vC in terms of u
syms t
syms u1(t)
syms u2(t)
syms u3(t)

xdot = cos(beta + theta)/cos(beta)*u1+CD*sin(theta)*u2; %H(1,:)*[u1;u2;0];
ydot = sin(beta + theta)/cos(beta)*u1-CD*cos(theta)*u2; %H(2,:)*[u1;u2;0];
thetadot = u2; %H(3,:)*[u1;u2;0];

vC = xdot*E1+ydot*E2;
% vD = vC+AD*thetadot*e2;

% calculating v1dot and v2dot
vCdot = diff(vC,t);


%% Finding the velocity of the mass center
syms u1dot u2dot real
% manually replacing the derivatives in the previous expressions by u1dot
% and u2dot
vCdot = [(cos(beta + theta)*u1dot)/cos(beta) + CD*sin(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) + CD*cos(theta)*u2*u2;
    (sin(beta + theta)*u1dot)/cos(beta) - CD*cos(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + CD*sin(theta)*u2*u2];

% manually squaring the previous expression
vcdot_squared = ((cos(beta + theta)*u1dot)/cos(beta) + CD*sin(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) + CD*cos(theta)*u2*u2)*((cos(beta + theta)*u1dot)/cos(beta) + CD*sin(theta)*u2dot - (sin(beta + theta)*u1*u2)/cos(beta) + CD*cos(theta)*u2*u2)+((sin(beta + theta)*u1dot)/cos(beta) - CD*cos(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + CD*sin(theta)*u2*u2)*((sin(beta + theta)*u1dot)/cos(beta) - CD*cos(theta)*u2dot + (cos(beta + theta)*u1*u2)/cos(beta) + CD*sin(theta)*u2*u2);
vcdot_squared = simplify(vcdot_squared);

%% Writing the Gibbs-Appell function
syms m real % the total mass of the system: m = m1+m2+m3
syms Izz real % moment of inertia of the vehicle about the z-axis passing through its mass center

S = m/2*vCdot_squared+Izz/2*u2dot*u2dot;

simplify(diff(S,u1dot))
simplify(diff(S,u2dot))

dSdu1dot = (m*(2*u1dot + 2*CD*cos(beta)^2*u2(t)^2 - CD*u2dot*sin(2*beta)))/(2*cos(beta)^2);
dSdu2dot = Izz*u2dot + CD^2*m*u2dot - CD*m*u1(t)*u2(t) - (CD*m*u1dot*sin(beta))/cos(beta);

A = [diff(dSdu1dot, u1dot), diff(dSdu1dot, u2dot);
    diff(dSdu2dot, u1dot), diff(dSdu2dot, u2dot)];

% A =
% 
% [              m/cos(beta)^2, -(CD*m*sin(2*beta))/(2*cos(beta)^2)]
% [-(CD*m*sin(beta))/cos(beta),                       m*CD^2 + Izz]


B = simplify([dSdu1dot;dSdu2dot]-A*[u1dot; u2dot]);
 
% B =
% 
%           (m*(2*u1dot + CD*cos(beta)^2*u2(t)^2 - CD*u2dot*sin(2*beta)))/cos(beta)^2
% 2*Izz*u2dot + 2*CD^2*m*u2dot - CD*m*u1(t)*u2(t) - (2*CD*m*u1dot*sin(beta))/cos(beta)


udot = A\B

simplify(udot)

% udot = 
% 
% (2*Izz*u1dot + CD^3*m*cos(beta)^2*u2(t)^2 + 2*CD^2*m*u1dot*cos(beta)^2 + Izz*CD*cos(beta)^2*u2(t)^2 - (CD^2*m*sin(2*beta)*u1(t)*u2(t))/2)/(m*cos(beta)^2*CD^2 + Izz)
%                                            (2*Izz*u2dot - CD*m*u1(t)*u2(t) + (CD^2*m*sin(2*beta)*u2(t)^2)/2 + 2*CD^2*m*u2dot*cos(beta)^2)/(m*cos(beta)^2*AD^2 + Izz)

