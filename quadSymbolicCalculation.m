%%
% clc
% clear all
% close all

phi = sym('phi','real');
theta = sym('theta','real');
psi = sym('psi','real');
eul = [phi; theta; psi];
Rquad = quadeul2rot(eul);
g = [0 0 sym('g','real')]';

% Velocity
vB = [sym('vBx','real'); sym('vBy','real'); sym('vBz','real')];
omega = [sym('p','real') sym('q','real') sym('r','real')]';
% R'*g
Tz = [0 0 sym('T_z','real')]';
% R*Uz
% quadSskew(omega)*vB
% cross(omega,vB)

% We compute the linearize model
RL = [1 -psi theta;psi 1 -phi;-theta phi 1]
TL = [1 0 theta;0 1 -phi;0 phi 1]
%%

% Rotation matrix with small theta and phi
R = Rotz(psi)*Roty(theta)*Rotx(phi);

% Simple rotation matrix
Rsimp = [ cos(psi), - sin(psi), theta*cos(psi) + phi*sin(psi);...
 sin(psi), cos(psi) , theta*sin(psi) - phi*cos(psi);...
   -theta,                           phi,                             1];

T = [1 phi*theta theta; 0 1 -phi;0 phi 1];

I = diag([sym('Ixx','real')*[1 1] sym('Izz','real')])
% Linear velocity
u = sym('u_b','real'); v = sym('v_b','real'); w = sym('w_b','real');
V_B = [ u v w]';
V_w = [sym('u_w','real') sym('v_w','real') sym('w_w','real')]';
lambda1 = sym('lambda_1','real');

vdot = Tz - lambda1*V_B - cross(omega,V_B) - Rsimp'*[0 0 g(3)]'-lambda1*Rsimp'*V_w
omegadot = -inv(I)*cross(I*omega,omega)
% pqrdot = 


function Rz = Rotz(phi)
Rz = [cos(phi) -sin(phi) 0;...
    sin(phi) cos(phi) 0;...
    0 0 1];
end

function Ry = Roty(theta)
Ry = [1 0 (theta);...
    0 1 0;...
    -(theta) 0 1];
end
function Rx = Rotx(psi)
Rx = [1 0 0;...
    0 1 -(psi);...
    0 (psi) 1];
end


