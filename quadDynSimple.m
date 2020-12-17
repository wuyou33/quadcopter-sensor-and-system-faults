function pddotOut = quadDynSimple(pdot,eulerdot,euler,U,V_W)
%quadDynSimple Simplified translation and rotation of quadcopter
%   euler: the Euler angles of quadcopter
%   U: The control signal around trim condition: U_z tau_phi tau_theta
%   tau_psi
%   V_W: wind velocity: u_w v_w w_w
param = quadparam;
u_w = V_W(1);v_w = V_W(2);w_w = V_W(3);
pddot = 1/param.m*([0 0 U(1)]' - param.lambda1*[pdot(1)+u_w; pdot(2)+v_w; pdot(3)+w_w] + ...
    param.g(3)*[euler(2); -euler(1); 0]);
eulddot = [U(2)/param.Jx; U(3)/param.Jy; U(3)/param.Jz]-...
    [param.alpha1*eulerdot(1); param.alpha1*eulerdot(2); param.alpha2*eulerdot(3)];
pddotOut = [pddot;eulddot];
end

