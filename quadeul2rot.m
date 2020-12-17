function R = quadeul2rot(euler)
%quadeul2rot convert euler angles 2 rotation
% Help
if nargin == 1
    opt = 'ZYX';
else
    if sum(strcmp(varargin{1},{'ZYZ';'ZYX'}))~=0
        opt = varargin{1};
    else
        opt = 'ZYX';
    end
end
phi = euler(1);theta = euler(2);psi = euler(3);
switch(opt)
    case 'ZYX'
        R = Rotz(psi)*Roty(theta)*Rotx(phi);
    case 'ZYZ'
        R = Rotz(phi)*Roty(theta)*Rotz(psi);
end
end

function Rz = Rotz(phi)
Rz = [cos(phi) -sin(phi) 0;...
    sin(phi) cos(phi) 0;...
    0 0 1];
end

function Ry = Roty(theta)
Ry = [cos(theta) 0 sin(theta);...
    0 1 0;...
    -sin(theta) 0 cos(theta)];
end
function Rx = Rotx(psi)
Rx = [1 0 0;...
    0 cos(psi) -sin(psi);...
    0 sin(psi) cos(psi)];
end


