%% Function to find the optimal path for the quadcopter.
clc

% The initial idea is to design the yaw angle (heading) as a increase and
% decrease ramp, which mean yaw rate is jump between r_max and -r_max, now
% we want to design for different yaw angle profile to figure out which one
% gives the best conditional number of the time-varying observability
% matrix

Ts = 0.005; T = 1; % We consider only 1 second, limited by GPS measurement
t = 0:Ts:T;
r_max = 0.1;
psi_max = 0.1;
psit = 0.2;
g = 9.81;
x(1) = r_max;
x(2) = psi_max;
rt = r_max;


g = 9.81;lambdam = 0.36/0.5;
S = [0 -1;1 0];
Rpsi = @(psi) [cos(psi) -sin(psi);sin(psi) cos(psi)];
Z2 = zeros(2,2);
Ar = @(psi, r) [-r*S Z2 Z2;g*S -r*S Z2;Z2 Rpsi(psi) Z2];
Cr_gps = @(psi, r) [Z2 Z2 eye(2);Z2 Rpsi(psi) Z2];

Aar = @(psi,r) [Ar(psi, r) [Z2 -eye(2);eye(2) Z2;Z2 Z2];zeros(4,6) zeros(4,4)];
Car_gps = @(psi,r) [Cr_gps(psi, r) [Z2 Z2;Z2 Z2];Z2 Z2 Z2 Z2 eye(2)];

As = @(psi,r) [-r*S Z2 Z2 Z2;g*S -lambdam*eye(2)-r*S Z2 -lambdam*Rpsi(psi)';...
    Z2 Rpsi(psi) Z2 Z2;Z2 Z2 Z2 Z2];
Cs_gps = @(psi, r) [Z2 Z2 eye(2) Z2;Z2 Rpsi(psi) Z2 Z2];
Cs_imu = @(psi,r) [Z2 lambdam*eye(2) Z2 lambdam*Rpsi(psi)'];

Aas = @(psi,r) [As(psi, r) [Z2 -eye(2); Z2 Z2;Z2 Z2;Z2 Z2];zeros(4,12)];
Cas_gps = @(psi,r) [Cs_gps(psi,r) zeros(4,4);zeros(2,10) eye(2)];
Cas_imu = @(psi,r) [Cs_imu(psi,r) eye(2) Z2];
% Testing the observability matrix
psifree = stateData(1:end,6);
psidotfree = stateData(1:end,18);
ind0 = 10;indbet = 1;indend = ind0+indbet*200;indend = ind0+2;
% for ii = 2:length(psifree)
%     indend = ind0+indbet*(ii-1);
%     [Oas,rankOas] = quadObsv(Aas, Cas_imu, Cas_gps, psifree(ind0:indbet:indend), psidotfree(ind0:indbet:indend));
%     condii(ii,:) = cond(Oas'*Oas);
%     ii
% end
[Oas,rankOas] = quadObsv(Aas, Cas_imu, Cas_gps, psifree(ind0:indbet:indend), psidotfree(ind0:indbet:indend));
for ii = ind0:indend
    Ct = [ Cas_imu(psifree(ii),psidotfree(ii));Cas_gps(psifree(ii),psidotfree(ii))];
    At = Aas(psifree(ii),psidotfree(ii));
    if ii == ind0
        productAtpre = eye(size(At,2));
        O = Ct*productAtpre;
        productAtpre = At*productAtpre;
    else
        [O,productAtpre] = quadObsvRecurive(At,Ct,O,productAtpre);
    end
    ii
end
ii = ind0;
C1 = [ Cas_imu(psifree(ii),psidotfree(ii));Cas_gps(psifree(ii),psidotfree(ii))];
A1 = Aas(psifree(ii),psidotfree(ii));
ii = ind0+1;
C2 = [ Cas_imu(psifree(ii),psidotfree(ii));Cas_gps(psifree(ii),psidotfree(ii))];
A2 = Aas(psifree(ii),psidotfree(ii));
ii = ind0+2;
C3 = [ Cas_imu(psifree(ii),psidotfree(ii));Cas_gps(psifree(ii),psidotfree(ii))];
A3 = Aas(psifree(ii),psidotfree(ii));
ii = ind0+3;
C4 = [ Cas_imu(psifree(ii),psidotfree(ii));Cas_gps(psifree(ii),psidotfree(ii))];
A4 = Aas(psifree(ii),psidotfree(ii));
O1 = C1;O2 = [C1;C2*A1];O3 = [C1;C2*A1;C3*A2*A1];O4 = [O3;C4*A3*A2*A1];

% [Oas2,rankOas2] = quadObsv(Aas, Cas_imu, Cas_gps, psinoise(ind0:indbet:indend), psidotnoise(ind0:indbet:indend));
% [Oar,rankOar] = quadObsv(Aar, Car_gps, [], psifree(ind0:indbet:indend), psidotfree(ind0:indbet:indend));
% 
% % cond(obsv(Aar(psifree(ind0),psidotfree(ind0)),Car_gps(psifree(ind0),psidotfree(ind0))))
% cond(obsv(Aas(psifree(ind0),psidotfree(ind0)),[Cas_gps(psifree(ind0),psidotfree(ind0));...
%     Cas_imu(psifree(ind0),psidotfree(ind0))]))
% cond(obsv(Aas(psinoise(ind0),psidotnoise(ind0)),[Cas_gps(psinoise(ind0),psidotnoise(ind0));...
%     Cas_imu(psinoise(ind0),psidotnoise(ind0))]))
% 
% rank(obsv(Aar(0,1),Car_gps(0,1)))
% rank(obsv(Aas(1,1),[Cas_gps(1,1);...
%     Cas_imu(1,1)]))

% [lambdaO, O, Aout, CN] = quadObsv([1 1 1],1);
% x0  =[0 0 0];
% options = optimoptions('fmincon'); 
% options.StepTolerance = 0.01;
% options.Display = 'none';
% for ii = 1:250
% ind0 = ii;
% X = fmincon(@(x) quadObsv2(x,ind0),x0,[],[],[],[],[0 0 0],[1 1 2*pi],[],options);
% Xm(ii,:) = X;
% fm(ii,:) = quadObsv(X,ind0);
% ii
% end
% quadObsv([1 1 1],ind0)





% function lambdaOut = quadObsv2(x,ind0)
% % Using psi = x(1)*cos(x(3)*t) + x(2)*sin(x(3)*t)
% % Then r = -x(1)*x(3)*sin(x(3)*t) + x(2)*x(3)*cos(x(3)*t)
% % Input: x is the variables of yaw angle profile
% % Nx, time instant when GPS and gyro bias available, can be varied
% % Aar: process model
% % Car1: 
% Ts = 0.005; T = 3; % We consider only 1 second, limited by GPS measurement
% t = 0:Ts:T;
% N = length(t);
% 
% if ind0>N-1/Ts
%     ind0 = 1;
% end
% psit =  x(1)*cos(x(3)*t) + x(2)*sin(x(3)*t);
% rt = -x(1)*x(3)*sin(x(3)*t) + x(2)*x(3)*cos(x(3)*t);
% g = 9.81;
% S = [0 -1;1 0];
% Rpsi = @(psi) [cos(psi) -sin(psi);sin(psi) cos(psi)];
% Z2 = zeros(2,2);
% Ar = @(psi, r) [-r*S Z2 Z2;g*S -r*S Z2;Z2 Rpsi(psi) Z2];
% Cr = @(psi, r) [Z2 Z2 eye(2);Z2 Rpsi(psi) Z2];
% Aar = @(psi,r) [Ar(psi, r) [Z2 -eye(2);eye(2) Z2;Z2 Z2];zeros(4,6) zeros(4,4)];
% Car = @(psi,r) [Cr(psi, r) [Z2 Z2;Z2 Z2];Z2 Z2 Z2 Z2 eye(2)];
% 
% C1 = Car(psit(1), rt(1));
% Aar1 = Aar(psit(1), rt(1));
% n = size(Aar1,1);
% m = size(C1,1);
% 
% Aout = eye(n);
% for ii = ind0:ind0+200
%     Aart = eye(n) + Ts*Aar(psit(ii), rt(ii));
%     Aout = Aart*Aout;
% end
% CN = Car(psit(N), rt(N));
% O = CN*Aout;
% lambdaOut = cond(obsv(Aout,CN));
% end




