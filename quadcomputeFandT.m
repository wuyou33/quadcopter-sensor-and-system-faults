function U = quadcomputeFandT(omega)
%quadcomputeTandF compute the forces and torques from motor speed
%   Detailed explanation goes here
param = quadparam;

U = param.T*omega.^2;
end

