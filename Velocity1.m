function [t3p, t4p] = Velocity1(A, t3, t4, t)
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

%********************************************
% finding equations for w2 and alpha2
w = pi()/180; %from degrees per second to radians per second (1 is given)
t2 = ((t2max - t2min)/2) * sin(w*t) + ((t2max + t2min)/2);
w2 = ((t2max - t2min)/2) * w * cos(w*t);
alp2 = -((t2max - t2min)/2) * w^2 * sin(w*t);

%Calculate the jacobian
J = [-r3 * sin(t3), -r4 * sin(t4); r3 * cos(t3), r4 * cos(t4)];

%solve for the first order KC values
KC1 = J\[r2*sin(t2) ; -r2*sin(t2)];

%solve for w3 and w4
t3p = KC1(1,1);
t4p = KC1(2,1);

end