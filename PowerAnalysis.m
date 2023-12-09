function torque = PowerAnalysis(A, thetamax, thetamin, w, t, r15, r5, r6, ...
    angle1, angle2, angle3, angle4, angle5, angle6, h6)

% Define link lengths(in meters)
A = A*10^(-3);
r1 = 4*A;
r2 = 5*A;
r3 = 2*A;
r4 = r2;
r5 = r5/1000;
r6 = r6/1000;
r26 = r3/2;
r15 = r15/1000;
h6 = h6/1000;
t1 = 0; % r1 is always horizontal to the right

% Define Kinematic Coefficients

J34 = [-r3*sin(angle3) -r4*sin(angle4);
    r3*cos(angle3) r4*cos(angle4)]; % Jacobian for first vector loop
sol341 = [r2*sin(angle2); 
    -r2*cos(angle2)]; % solution vector for first kinematic coefficients

theta341 = J34\sol341; % solve for first order kinematic coefficients for
                       % theta 3, theta 4

theta31 = theta341(1);
theta41 = theta341(2);

sol342 = [r2*cos(angle2) + theta31^2*r3*cos(angle3) + theta41^2*r4*cos(angle4);
    r2*sin(angle2) + theta31^2*r3*sin(angle3) + theta41^2*r4*sin(angle4)]; % solution vector for the second order kinematic coefficients

theta342 = J34\sol342; % solve for the second order kiinematic coefficients for theta 3, theta 4

theta32 = theta342(1);
theta42 = theta342(2);


% second loop

J56 = [-r5*sin(angle5) -r6*sin(angle6);
    r5*cos(angle5) r6*cos(angle6)]; % Jacobian for the second loop 
sol561 = [r2*sin(angle2) + r26*theta31*sin(angle3) + r15*theta41*sin(angle4); 
    -r2*cos(angle2) - r26*theta31*cos(angle3) - r15*theta41*cos(angle4)]; % solution vector 

theta561 = J56\sol561; % solve for first order kinematic coefficients for
                       % theta 5, theta 6

theta51 = theta561(1);
theta61 = theta561(2);

sol562 = [(r2*cos(angle2) + r26*theta32*sin(angle3) + r26*theta31^2*cos(angle3) + ...
    r6*theta61^2*cos(angle6) + r5*theta51^2*cos(angle5) + r15*theta42*sin(angle4) + ...
    r15*theta41^2*cos(angle4));
    (r2*sin(angle2) - r26*theta32*cos(angle3) + r26*theta31^2*sin(angle3) + ...
    r6*theta61^2*sin(angle6) + r5*theta51^2*sin(angle5) - r15*theta42*cos(angle4) + ...
    r15*theta41^2*sin(angle4));];

theta562 = J56\sol562; % solve for second order kinematic coefficients for
                       % theta 5, theta 6

theta52 = theta562(1);
theta62 = theta562(2);


% Define inputs
% Calculate mass
rho = 2698.9; % kg/m^3 (density)
% rho = rho / 1000^3; % kg/mm^3
dL = 30/1000; % thickness of link
wL = 30/1000; % width of link

% Actual lengths
r1a = r1 + wL;
r2a = r2 + wL;
r3a = r3 + wL;
r4a = r4 + wL;
r5a = r5 + wL;
r6a = r6 + wL;

% Mass in kg
m1 = rho*r1a*wL*dL;
m2 = rho*r2a*wL*dL;
m3 = rho*r3a*wL*dL;
m4 = rho*r4a*wL*dL;
m5 = rho*r5a*wL*dL;
m6 = rho*r6a*wL*dL + 2*rho*(h6-wL)*wL*dL;

% m1 = rho*(r1+wL)*wL*dL; % kg
% m2 = rho*(r2+wL)*wL*dL;
% m3 = rho*(r3+wL)*wL*dL;
% m4 = rho*(r4+wL)*wL*dL;
% m5 = rho*(r5+wL)*wL*dL;
% m6 = rho*(r6+wL)*wL*dL;

% Calculate moments of inertia
i1 = 0;
i2 = (1/12)*m2*((wL^2)+(r2)^2);
i3 = (1/12)*m3*((wL^2)+(r3)^2);
i4 = (1/12)*m4*((wL^2)+(r4)^2);
i5 = (1/12)*m5*((wL^2)+(r5)^2);
i6 = (1/12)*m6*((wL^2)+(r6)^2);

%% Repeat for 3 periods of oscillation

% Position Analysis of Vector Loop 1

% Position Analysis of Vector Loop 2
 
% Velocity Analysis of Vector Loop 1

% Velocity Analysis of Vector Loop 2

% Acceleration Analysis of Vector Loop 1

% Acceleration Analysis of Vector Loop 2

% Find KC of center of mass of link 2, 3, 4, 5, 6

% Naming convention - xCoMr2 is the x coords of the center of mass of link 2
%                     x1CoMr2 is the x rate of change of the center of mass of link 2
%                     x2CoMr2 is the x acceleration of the center of mass of link 2
%                     CoMa2 is the angle of the vector which points to the center of mass of link 2
%                     All of the rest of the CoM terms are named in this convention for the different links

xCoMr2 = 0.5*r2*cos(angle2);
yCoMr2 = 0.5*r2*sin(angle2);
x1CoMr2 = -0.5*r2*sin(angle2);
y1CoMr2 = 0.5*r2*cos(angle2);
x2CoMr2 = -xCoMr2;
y2CoMr2 = -yCoMr2;
CoMa2 = angle2;

xCoMr3 = (r2*cos(angle2))+(0.5*r3*cos(angle3));
yCoMr3 = (r2*sin(angle2))+(0.5*r3*sin(angle3));
x1CoMr3 = (-r2*sin(angle2))-(0.5*theta31*r3*sin(angle3));
y1CoMr3 = (r2*cos(angle2))-(0.5*theta31*r3*cos(angle3));
x2CoMr3 = (-r2*cos(angle2))-(0.5*theta32*r3*sin(angle3))-(0.5*(theta31)^2*r3*cos(angle3));
y2CoMr3 = (-r2*sin(angle2))+(0.5*theta32*r3*cos(angle3))-(0.5*(theta31)^2*r3*sin(angle3));
CoMa3 = atan(yCoMr3/xCoMr3);


xCoMr4 = (r1*cos(angle1))-(0.5*r4*cos(angle4));
yCoMr4 = (r1*sin(angle1))+(0.5*r4*sin(angle4));
x1CoMr4 = (0.5*r4*theta41*sin(angle4));
y1CoMr4 = (-0.5*r4*theta41*cos(angle4));
x2CoMr4 = (0.5*theta42*r4*sin(angle4))+(0.5*(theta41)^2*r4*cos(angle4));
y2CoMr4 = (-0.5*theta42*r4*cos(angle4))+(0.5*(theta41)^2*r4*sin(angle4));
CoMa4 = atan((yCoMr4)/(xCoMr4));

xCoMr5 = (r1*cos(angle1))-(r15*cos(angle4))-(0.5*r5*cos(angle5));
yCoMr5 = (r1*sin(angle1))-(r15*sin(angle4))-(0.5*r5*sin(angle5));
x1CoMr5 = (r15*theta41*sin(angle4))+(0.5*r5*theta51*sin(angle5));
y1CoMr5 = (-r15*theta41*cos(angle4))-(0.5*r5*theta51*cos(angle5));
x2CoMr5 = (r15*theta42*sin(angle4))+(r15*(theta41)^2*cos(angle4))+(0.5*r5*theta52*sin(angle5))+(0.5*r5*(theta51)^2*cos(angle5));
y2CoMr5 = (-r15*theta42*cos(angle4))+(r15*(theta41)^2*sin(angle4))-(0.5*r5*theta52*cos(angle5))+(0.5*r5*(theta51)^2*sin(angle5));
CoMa5 = atan(yCoMr5/xCoMr5);


xCoMr6 = (r1*cos(angle1))-(r15*cos(angle4))-(r5*cos(angle5))-(0.5*r6*cos(angle6))+(r6*cos(angle6+pi));
yCoMr6 = (r1*sin(angle1))-(r15*sin(angle4))-(r5*sin(angle5))-(0.5*r6*sin(angle6))+(r6*sin(angle6+pi));
x1CoMr6 = (r15*theta41*sin(angle4))+(r5*theta51*sin(angle5))+(0.5*r6*theta61*sin(angle6))-(r6*theta61*sin(angle6+pi));
y1CoMr6 = (-r15*theta41*cos(angle4))-(r5*theta51*cos(angle5))-(0.5*r6*theta61*cos(angle6))+(r6*theta61*cos(angle6+pi));
x2CoMr6 = (r15*theta42*sin(angle4))+(r15*(theta41)^2*cos(angle4))+(r5*theta52*sin(angle5))+(r5*(theta51)^2*cos(angle5))+(0.5*r6*theta62*sin(angle6))+(0.5*r6*(theta61)^2*cos(angle6))-(r6*theta62*sin(angle6+pi))-(r6*(theta61)^2*cos(angle6+pi));
y2CoMr6 = (-r15*theta42*cos(angle4))+(r15*(theta41)^2*sin(angle4))-(r5*theta52*cos(angle5))+(r5*(theta51)^2*sin(angle5))-(0.5*r6*theta62*cos(angle6))+(0.5*r6*(theta61)^2*sin(angle6))+(r6*theta62*cos(angle6+pi))-(r6*(theta61)^2*sin(angle6+pi));
CoMa6 = atan(yCoMr6/xCoMr6);


%  Given constant omega
% w = 1 * pi / 180;
thetamax = thetamax*pi/180;
thetamin = thetamin*pi/180;

%  Omega and Alpha arrays (link 2)
omega2 = w / 2 * (thetamax - thetamin) * cos(w * t);
alpha2 = -w^2 / 2 * (thetamax - thetamin) * sin(w * t);

% Calculate the Torque Needed at any Time 
% P = T2 * w2 = sum(Aj * w2 * alpha2) + sum(Bj * w2^2) + sum(mj * g * y1CoMrj * w2)
% Can cancel w2 on both sides
% T2 = sum(Aj * alpha2) + sum(Bj * w2) + sum(mj * g * y1CoMrj)

% Calculate the Change in Kinetic Energy dT/dt using the Aj and Bj functions
A2 = Aj(m2,x1CoMr2,y1CoMr2, i2, 1);
A3 = Aj(m3,x1CoMr3,y1CoMr3, i3, theta31);
A4 = Aj(m4,x1CoMr4,y1CoMr4, i4, theta41);
A5 = Aj(m5,x1CoMr5,y1CoMr5, i5, theta51);
A6 = Aj(m6,x1CoMr6,y1CoMr6, i6, theta61);

B2 = Bj(m2, x1CoMr2, x2CoMr2, y1CoMr2, y2CoMr2, i2, 1, 0);
B3 = Bj(m3, x1CoMr3, x2CoMr3, y1CoMr3, y2CoMr3, i3, theta31, theta32);
B4 = Bj(m4, x1CoMr4, x2CoMr4, y1CoMr4, y2CoMr4, i4, theta41, theta42);
B5 = Bj(m5, x1CoMr5, x2CoMr5, y1CoMr5, y2CoMr5, i5, theta51, theta52);
B6 = Bj(m6, x1CoMr6, x2CoMr6, y1CoMr6, y2CoMr6, i6, theta61, theta62);

A = (A2 + A3 + A4 + A5 + A6) * omega2 * alpha2;
B = (B2 + B3 + B4 + B5 + B6) * omega2 ^ 3;

dTdt = A + B;

% Calculate the Change in Gravitational Potential Energy dU/dt

Ugrav2 = m2 * 9.81 * y1CoMr2;
Ugrav3 = m3 * 9.81 * y1CoMr3;
Ugrav4 = m4 * 9.81 * y1CoMr4;
Ugrav5 = m5 * 9.81 * y1CoMr5;
Ugrav6 = m6 * 9.81 * y1CoMr6;

dUdt= (Ugrav2 + Ugrav3 + Ugrav4 + Ugrav5 + Ugrav6) * omega2;

% Calculate Torque at a Time t

torque = dUdt + dTdt;

end



