function omega = quadcomputeMotorSpeed(U)
%quadcomputeMotorSpeed Compute motor speed from required torques and forces

%   Detailed explanation goes here
param = quadparam;
% Control signal
omega2 = inv(param.T)*U;
omega = sqrt(omega2);
end

