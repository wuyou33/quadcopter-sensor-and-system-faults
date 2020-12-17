function param = quadparam
%quadparam Initialize the quad parameters
%   Detailed explanation goes here
param.m = 0.5;
param.g = [0 0 9.81]';
param.cd = 0.32;
param.cf = 1e-4;
param.L = 0.5;
param.ct = 1e-5;
param.Jx = 0.002;
param.Jy = 0.002;
param.Jz = 0.003;

param.lambda1 = 0.36;
param.lambda2 = 0.6;

param.alpha1 = 0.002;
param.alpha2 = 0.005;
param.alphaPhi = 0.000;

param.dt = 5*1e-3;
param.tfinal = 120;

param.Jmatrix = diag([param.Jx param.Jy param.Jz]);
param.invJmatrix = param.Jmatrix\eye(3);

param.dragLin = diag([param.lambda1 param.lambda1 param.lambda2]);
param.dragRot = diag([param.alpha1 param.alpha1 param.alpha2]);

param.T = [param.cf*[1 1 1 1];...
    param.L*param.cf*[-1 1 1 -1];...
    param.L*param.cf*[-1 -1 1 1];...
    param.ct*[-1 1 -1 1]];

param.invT = param.T\eye(4);
end

