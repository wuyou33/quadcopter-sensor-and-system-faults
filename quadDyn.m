function pddotOut = quadDyn(U,pe,pedot,eul,pqr)
%quadDyn main dynamic model
% [pdot,euldot,pqrdot] = quadDyn(px,pxdot,eul,euldot,pqr,param)
% Define the state: position pe(x,y,z) and euler (phi,theta,psi)
% Define the derivative of state: pedot, eulerdot xdot

% Initial the parameters of quadcopter
param = quadparam;

% Control signal
% U = [0 0 0 param.m*param.g(3)];

Uz = U(1); 
tau_phi = U(2);
tau_theta = U(3);
tau_psi = U(4);


% We have the euler dot, need to convert to pqr
p = pqr(1);q = pqr(2); r = pqr(3);

% The dynamical of pqr
Jx = param.Jx;
Jy = param.Jy;
Jz = param.Jz;

pqrdot = [(Jy-Jz)/Jx*q*r; (Jz-Jx)/Jy*p*r; (Jx-Jy)/Jz*p*q] +...
    [tau_phi/Jx; tau_theta/Jy;tau_psi/Jz];

% Compute the body velocity
R = quadeul2rot(eul);
vB = R'*pedot;
pddot = 1/param.m*( [0 0 Uz]' - diag([param.cd param.cd param.cd])*vB ) - R'*param.g;

% % Sensor measurement
% acc_meas = 1/param.m*( [0 0 -Uz]' - diag([param.cd param.cd param.cd])*vB );
% pddot = acc_meas;
% % Translation model
% vBdot =  acc_meas + param.g*[-sin(theta); cos(theta)*sin(phi); cos(phi)*cos(theta)];
pddotOut = [pddot; pqrdot];
end

