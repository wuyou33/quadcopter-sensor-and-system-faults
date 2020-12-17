function Rpqr = quadRpqr2euldot(euler)
%quadRqr2euldot compute the 3x3 relation, from pqr to euler dot
%   If we want to call the function from SIMULINK, the output should be a
%   vector

phi  = euler(1);
theta = euler(2);
Rpqr = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
    0 cos(phi) -sin(phi);...
    0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
end