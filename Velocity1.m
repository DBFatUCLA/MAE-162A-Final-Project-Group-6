function [t3p, t4p, t3pp, t4pp] = Loop1(A, t3, t4, t)
% Velocity analysis of vector loop one

% A = 200;
r1 = 4*A; %mm
r2 = 5*A;
r3 = 2*A;
r4 = r2;
% t1 = 0; %rad
% 
% % Set rightmost deflection angles
% t2 = atan(3/4);
% t3 = pi()/2;
% t4 = 3*pi()/2;
% r26 = r3/2;
% t15 = t4;
% t26 = t3;
% r5 = 850;

%Calculate the jacobian
J = [-r3 * sin(t3), -r4 * sin(t4); r3 * cos(t3), r4 * cos(t4)];

%********************************************
% finding equations for w2 and alpha2
w = pi()/180; %from degrees per second to radians per second (1 is given)

alp2 = -((t2max - t2min)/2) * w^2 * sin(w*t);

%solve for the first order KC values
KCv = J\[r2*sin(t2) ; -r2*sin(t2)];

% isolate first order kinematic coefficients
t3p = KCv(1,1);
t4p = KCv(2,1);

%% solve for acceleration
KCa = J \ [r2*cos(t2) + t3p^2 * r3*cos(t3) + t4p^2 * r4*cos(t4) ; 
           r2*sin(t2) + t3p^2 * r3 *sin(t3) + t4p^2 * r4 * sin(t4)];

% isolate second order kinematic coefficients
t3pp = KCa(1,1);
t4pp = KCa(2,1);

end