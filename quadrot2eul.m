function euler = quadrot2eul(R,varargin)
%quadrot2eul Function to find euler angles from rotation matrix
% Two options: can be 0pi or 0-pi
% We have to check R(3,3) ~= 1
if R(3,3) ~=1 % theta~=pi/2
    if nargin == 1
        option = '0pi';
    else
        option = varargin{1};
    end
    switch(option)
        case '0pi'
            phi = atan2( R(2,3), R(1,3));
            theta = atan2( sqrt( R(1,3)^2 + R(2,3)^2), R(3,3));
            psi = atan2( R(3,2), -R(3,1));
        case '0-pi'
            phi = atan2( -R(2,3), -R(1,3));
            theta = atan2( -sqrt( R(1,3)^2 + R(2,3)^2), R(3,3));
            psi = atan2( -R(3,2), R(3,1));
    end
    euler = [phi, theta, psi]';
else
    thetapsi = atan2(R(2,1),R(2,2));
    psi = 0;
    phi = thetapsi-psi;
    theta = acos(R(3,3));
    euler = [phi theta psi]';
end

