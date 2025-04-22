% obtaining the Gibbs-Appell equations of motion

%% Defining the quasi-coordinates
% finding the matrix H such that
% [xdot; ydot; thetadot] = H*[u1; u2; u3]

syms l gamma real
syms theta(t)
syms xdot ydot thetadot real

E1 = [1,0];
E2 = [0,1];

e1 = cos(theta)*E1+sin(theta)*E2;
e2 = -sin(theta)*E1+cos(theta)*E2;

vc = xdot*E1+ydot*E2;
v3 = vc+l*thetadot*e2;

u1 = xdot*cos(theta)+ydot*sin(theta); %dot(vc,e1);
u2 = thetadot;
u3 = l*thetadot - tan(gamma)*(xdot*cos(theta) + ydot*sin(theta)) + ydot*cos(theta) - xdot*sin(theta); %dot(vc,e2)+l*thetadot-tan(gamma)*dot(vc,e1);

J = [diff(u1,xdot), diff(u1, ydot), diff(u1,thetadot);
    diff(u2,xdot), diff(u2, ydot), diff(u2,thetadot);
    diff(u3,xdot), diff(u3, ydot), diff(u3,thetadot)];

H = simplify(inv(J));


% writing v1 and v2 in terms of u
syms t
syms u1(t)
syms u2(t)
syms u3(t)

xdot = cos(gamma + theta)/cos(gamma)*u1+l*sin(theta)*u2; %H(1,:)*[u1;u2;0];
ydot = sin(gamma + theta)/cos(gamma)*u1-l*cos(theta)*u2; %H(2,:)*[u1;u2;0];
thetadot = u2; %H(3,:)*[u1;u2;0];

vc = xdot*E1+ydot*E2;
% v3 = vc+l*thetadot*e2;

% calculating v1dot and v2dot
vcdot = diff(vc,t);


%% Finding the velocity of the mass center
syms u1dot u2dot real
% manually replacing the derivatives in the previous expressions by u1dot
% and u2dot
vcdot = [(cos(gamma + theta)*u1dot)/cos(gamma) + l*sin(theta)*u2dot - (sin(gamma + theta)*u1*u2)/cos(gamma) + l*cos(theta)*u2*u2;
    (sin(gamma + theta)*u1dot)/cos(gamma) - l*cos(theta)*u2dot + (cos(gamma + theta)*u1*u2)/cos(gamma) + l*sin(theta)*u2*u2];

% manually squaring the previous expression
vcdot_squared = ((cos(gamma + theta)*u1dot)/cos(gamma) + l*sin(theta)*u2dot - (sin(gamma + theta)*u1*u2)/cos(gamma) + l*cos(theta)*u2*u2)*((cos(gamma + theta)*u1dot)/cos(gamma) + l*sin(theta)*u2dot - (sin(gamma + theta)*u1*u2)/cos(gamma) + l*cos(theta)*u2*u2)+((sin(gamma + theta)*u1dot)/cos(gamma) - l*cos(theta)*u2dot + (cos(gamma + theta)*u1*u2)/cos(gamma) + l*sin(theta)*u2*u2)*((sin(gamma + theta)*u1dot)/cos(gamma) - l*cos(theta)*u2dot + (cos(gamma + theta)*u1*u2)/cos(gamma) + l*sin(theta)*u2*u2);
vcdot_squared = simplify(vcdot_squared);

%% Writing the Gibbs-Appell function
syms m real % the total mass of the system: m = m1+m2+m3
syms Izz real % moment of inertia of the vehicle about the z-axis passing through its mass center

S = m/2*vcdot_squared+Izz/2*u2dot*u2dot;

simplify(diff(S,u1dot))
simplify(diff(S,u2dot))

dSdu1dot = (m*(2*u1dot + 2*l*cos(gamma)^2*u2(t)^2 - l*u2dot*sin(2*gamma)))/(2*cos(gamma)^2);
dSdu2dot = Izz*u2dot + l^2*m*u2dot - l*m*u1(t)*u2(t) - (l*m*u1dot*sin(gamma))/cos(gamma);

A = [diff(dSdu1dot, u1dot), diff(dSdu1dot, u2dot);
    diff(dSdu2dot, u1dot), diff(dSdu2dot, u2dot)];

% A =
% 
% [              m/cos(gamma)^2, -(l*m*sin(2*gamma))/(2*cos(gamma)^2)]
% [-(l*m*sin(gamma))/cos(gamma),                         m*l^2 + I_zz]


B = simplify([dSdu1dot;dSdu2dot]-A*[u1dot; u2dot]);
 
% B =
% 
%           (m*(2*u1dot + l*cos(gamma)^2*u2(t)^2 - l*u2dot*sin(2*gamma)))/cos(gamma)^2
% 2*I_zz*u2dot + 2*l^2*m*u2dot - l*m*u1(t)*u2(t) - (2*l*m*u1dot*sin(gamma))/cos(gamma)


udot = A\B

simplify(udot)

% udot = 
% 
% (2*I_zz*u1dot + l^3*m*cos(gamma)^2*u2(t)^2 + 2*l^2*m*u1dot*cos(gamma)^2 + I_zz*l*cos(gamma)^2*u2(t)^2 - (l^2*m*sin(2*gamma)*u1(t)*u2(t))/2)/(m*cos(gamma)^2*l^2 + I_zz)
%                                              (2*I_zz*u2dot - l*m*u1(t)*u2(t) + (l^2*m*sin(2*gamma)*u2(t)^2)/2 + 2*l^2*m*u2dot*cos(gamma)^2)/(m*cos(gamma)^2*l^2 + I_zz)

