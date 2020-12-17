%% CHECK THE MODEL
clc
clear variables
close all

omega = rand(3,1);
eul = 0.1*rand(3,1);
pe0 = [0 0 0]';
p = omega(1);q = omega(2);r = omega(3);
param = quadparam;
I = diag([param.Jx param.Jy param.Jz]);

% Same
Jx = param.Jx;Jy = param.Jy;Jz = param.Jz;
-quadSskew(omega)*(I*omega);
I*[(Jy-Jz)/Jx*q*r; (Jz-Jx)/Jy*p*r; (Jx-Jy)/Jz*p*q];

% Now we test the translation model

% 01. compute the rotation matrix
vB = rand(3,1);
R = quadeul2rot(eul);
g = [0 0 9.81]';
Uz = [0 0 rand]';
vBdot = Uz - R'*g - quadSskew(omega)*vB;

pdot = R*vB;
% pddot = quadSskew(omega)*R*vB
pddot = R*Uz - g;
Rdot = R*quadSskew(omega);

vdot = R'*(pddot - Rdot*vB);

% dt = 1e-10;
% Rpqr = quadpqr2euldot(eul);
% euldot = Rpqr*omega;
% dt = 1e-6;
% euln = eul+euldot*dt;
% Rn = quadeul2rot(euln);
% Rdotn = (Rn-R)/dt;

T = 1e-3;tfinal = 10;
dt = 1e-8;
euln = eul+dt*randn(3,1);
euldot = (euln-eul)/dt;
Rn = quadeul2rot(euln);
Rdotn = (Rn-R)/dt;
Rpqr = quadRpqr2euldot(eul);
Somega = R'*Rdotn;
omegan = [Somega(3,2) -Somega(3,1) Somega(2,1)];
pqr = Rpqr\euldot;

% Correct pqr?
R'*Rdotn-quadSskew(pqr);
pddot = quadDyn([param.m*param.g(3) 0 0 0],[0 0 0]',[0 0 0]',[0 0 0]',[0 0 0]');
disp('Done')

% Control tuning
Gphi = tf(1/param.Jx,[1 param.alpha1/param.Jx 0]);
Gtheta = tf(1/param.Jy,[1 param.alpha1/param.Jy 0]);
Gpsi = tf(1/param.Jz,[1 param.alpha2/param.Jz 0]);

% Tune for the orientation
% eulddot = [U(2)/param.Jx; U(3)/param.Jy; U(3)/param.Jz]-...
%     [param.alpha1*eulerdot(1)/param.Jx; param.alpha1*eulerdot(2)/param.Jy; param.alpha2*eulerdot(3)/param.Jz];
% ameas= 1/param.m*([0 0 U(1)]' - param.lambda1*[pdot(1)+u_w; pdot(2)+v_w;...
%     param.lambda2/param.lambda1*(pdot(3)+w_w)]);

Uzsat = 100*3*param.m*param.g(3);
Uphisat = 100*param.m*param.g(3)*0.5;
Uthetasat = 100*param.m*param.g(3)*0.5;
Upsisat = 100*param.m*param.g(3)*0.1;
s = tf('s');
Gphi = tf([0 0 1/param.Jx],[1 param.alpha1/param.Jx 0]);

Gx = tf(1,[1 param.lambda1/param.m 0]);
Gz = tf(1,[1 param.lambda2/param.m 0]);

% PID controller for position
param = quadparam;
% Kpx = 1.363;
% Kix = 0.44189;
% Kdx = 4.1557;
% Nx = 1/0.002099;

Kpx = 1.2418;
Kix = 0.16881;
Kdx = 1.7904;
Nx = 1/0.002099;


Kpy = 1.2418;
Kiy = 0.16881;
Kdy = 1.7904;
Ny = 1/0.002099;


Kpz = 1.2418;
Kiz = 0.16881;
Kdz = 1.7904;
Nz = 1/0.002099;

% PID control for the orientation
Kpphi = 0.0439*1;
Kiphi = 0.0098*0;
Kdphi = 0.0478*1;
Tfphi = 0.0061486*1;
Nphi = 1/Tfphi;

Kptheta = 0.0439*1;
Kitheta = 0.0098*0;
Kdtheta = 0.0478*1;
Tftheta = 0.0061486;
Ntheta = 1/Tftheta;

% raising 0.153s
Kppsi = 0.11038*1;
Kipsi = 0.045028*.0;
Kdpsi = 0.067644*1;
Tfpsi = 0.0062821;
Npsi = 1/Tfpsi;


% % We test simple tuning
Kpx = 5;
Kix = 1;
Kdx = 15;


Kpy = 5;
Kiy = 1;
Kdy = 15;



Kpz = 3.5;
Kiz = 1.5;
Kdz = 5.5;

% % PID control for the orientation
% Kpphi = 1;
% Kiphi = 0.;
% Kdphi = 1;
%
%
% Kptheta = 1;
% Kitheta = 0.0;
% Kdtheta = 1;
%
%
% % raising 0.153s
% Kppsi = 1;
% Kipsi = 0.0;
% Kdpsi = 1;



%

% SIMULATION COEFFICIENT
dt = 0.005; % 200Hz
tfinal = 120;
clock = 0:dt:tfinal;
N = length(clock);
% Design the velocity
Tz = [0 2 2 4 4 38 38 40 40 42]*2;
Vz = [0 0 2 2 0  0  0  0  0  0]*0;

% The square trajectory
Tx = [0 2 4 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 120];
Vx = [0 0 0  0  1  1  0  0  0  0  0  0  0  0 -1 -1  0  0  0  0  0   0   0   0   0];
Vy = [0 0 0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0 -1 -1   0   0   0   0];
Tpsi = [0 25 25 50 50 75 75 120];
psiref = [0  0 pi/2 pi/2 pi pi 3*pi/2 3*pi/2];
Tpsi = clock;

% Ramp+constant
ccpsi = pi/9*[1/500*(1:500) ones(1,1500) ];
psiref = [0 zeros(1,12000) ccpsi pi/9+ccpsi 2*pi/9+ccpsi 3*pi/9+ccpsi 4*pi/9+ccpsi,...
    5*pi/9+ccpsi];

% Instant jump
psiref = pi/9*[zeros(1,12001) ones(1,4000) 1.5*ones(1,4000) 1.5+0.5/2000*(1:2000) 2*ones(1,2000)];

% Ramp back and forth
ccpsi = pi/9*[1/300*(1:300) ones(1,400)  1-1/300*(1:300) ];
ccpsi = pi/5*[1/300*(1:300)  1-1/300*(1:300) ];

psiref = [0 kron(ones(1,20),[ccpsi -ccpsi])];

% % psiddot is rectangular
% ccii = 0.1*[ones(1,500) -ones(1,500) -ones(1,500) ones(1,500)];
% psiddot = [0 kron(ones(1,12),ccii)];
% psidot = filter(param.dt, [1 -1], psiddot);
% psiref = filter(param.dt, [1 -1], psidot);

% Sinusoidal
% ccpi = pi/9*(sin(2*pi*3.6/200*(1:length(clock))) + sin(2*pi*4.6/200*(1:length(clock))) + ...
%     sin(2*pi*5.6/200*(1:length(clock))));
% psiref = ccpi;

% psiref = pi/9*sin(2*pi/5*param.dt*(1:length(clock))).*cos(2*pi/4*param.dt*(1:length(clock))) + ...
%     pi/9*sin(2*pi/3*param.dt*(1:length(clock))).*cos(2*pi/2*param.dt*(1:length(clock))) + ...
%     pi/9*sin(2*pi/2*param.dt*(1:length(clock))).*cos(2*pi/1.5*param.dt*(1:length(clock)));

% The line trajectory
% TxL = [0 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
% VxL = [0 00 01 01 00 -1 -1 00 01 01 00 -1 -1 00 01 01 00 -1 -1  00  01  01  00   0];
% VyL = [0 00 01 01 00 -1 -1 00 01 01 00 -1 -1 00 01 01 00 -1 -1  00  01  01  00   0];

% % Line with shape change in velocity
% TxL = [0 10 10 25 25 40 40 55 55 70 70 85 85 100 100 115 115 120];
% VxL = [0 00 01 01 -1 -1 01 01 -1 -1 01 01 -1  -1  01  01  00  00];
% VyL = [0 00 01 01 -1 -1 01 01 -1 -1 01 01 -1  -1  01  01  00  00];

% % Line with constant veloctiy
% TxL = [0 10 10 120];
% VxL = [0 0 1 1];
% VyL = VxL;

% Sin velocity
GxL = tf(1,[1 1]);
GdxL = c2d(GxL,param.dt,'zoh');
[bxL,axL] = tfdata(GdxL,'v');
vfinal = 0;tmidd = 0;
TxL1 = [0 tmidd]; TxL2 = tmidd:param.dt:param.tfinal;
VxL1 = [0 vfinal];
gainV = 0.2;
VxL2 = vfinal + gainV*sin(2*pi/5*(TxL2-tmidd)).*cos(2*pi/4*(TxL2-tmidd)) + ...
    gainV*sin(2*pi/3*(TxL2-tmidd)).*cos(2*pi/2*(TxL2-tmidd)) + ...
    gainV*sin(2*pi/2*(TxL2-tmidd)).*cos(2*pi/1.5*(TxL2-tmidd));

VyL1 = [0 vfinal];
VyL2 = VxL2;


% VxL2 = VyL2;
TxL = [TxL1 TxL2];
VxL = [VxL1 VxL2];
VyL = [VyL1 VyL2];

TpsiL =   [0 25   25   40 40 55     55     70   70   85    85    100  100 115   115  120];
psirefL = [0  0 pi/2 pi/2 pi pi 3*pi/2 3*pi/2 2*pi 2*pi 3*pi/2 3*pi/2  pi  pi  pi/2 pi/2];
TxQ = [0 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
VxQ = [0  0  1  1  0  0  0  0 -1 -1  0  0  0  0  1  1  0  0  0   0   0   0   0   0];
VyQ = [0  0  0  0  0  1  1  0  0  0  0 -1 -1  0  0  0  0  1  1   0   0   0   0   0];
VxQsim = [TxQ' VxQ'];
VyQsim = [TxQ' VyQ'];


Teffpsi = [0 30 30 40 40 80 80 120];
effz = ones(size(Tpsi));
% effz = [1 1 1 1 0.8 0.8 0.8 0.8];
effPhi = ones(size(Tpsi));
effTheta = ones(size(Tpsi));
% effTheta = [1 1 1 1 0.8 0.8 0.8 0.8];
effPsi = ones(size(Tpsi));

effz = [Tpsi' effz'];
effPhi = [Tpsi' effPhi'];
effTheta = [Tpsi' effTheta'];
effPsi = [Tpsi' effPsi'];


Tphi =   [0 6  10  30  80  100 120];
Phidot = [0 0   0 0.1  0.1   0.1   0.1];

Phidotsim = [Tphi' Phidot'];

% line or square trajectory
pathLine = 1;
if pathLine == 1
    Vxsim = [TxL' VxL'];
    Vysim = [TxL' VyL'];
    Vzsim = [Tz' Vz'];
    psirefIn = [Tpsi' psiref'];
else
    Vxsim = [Tx' Vx'];
    Vysim = [Tx' Vy'];
    Vzsim = [Tz' Vz'];
    psirefIn = [Tpsi' psiref'];
    
end


% We have some choice
% Path selection 1 square 2 circle 3 8 shape
pathSelect = 1;
gainPsi = 1;
gainWind = 1;
gainAcc = 1;
gainGyro = 1;




% SIMULATION OF WIND TURBALENCE MODEL
sigmau = 1.1006;sigmav = 0.6576;sigmaw = 0.4168;
h = 8;
Lu = 280*(h/500)^.35;
Lv = 140*(h/500)^.48;
Lw = 0.35*h;
Va = sqrt(2);
Hu = tf(sigmau*sqrt(2*Va/Lu),[1 Va/Lu]);
Hv = tf(sigmav*sqrt(3*Va/Lv)*[1 Va/sqrt(3)/Lv],[1 2*Va/Lv Va^2/Lv^2]);
% Hv = Hu;
Hw = tf(sigmaw*sqrt(3*Va/Lw)*[1 Va/sqrt(3)/Lw],[1 2*Va/Lw Va^2/Lw^2]);
[bHu,aHu] = tfdata(Hu,'v');
[bHv,aHv] = tfdata(Hv,'v');
[bHw,aHw] = tfdata(Hw,'v');


windLevel = 1;
turuw = [clock' sqrt(windLevel)*randn(N,1)];
turvw = [clock' sqrt(windLevel)*randn(N,1)];
turww = [clock' sqrt(windLevel)*randn(N,1)];


% slow wind in x, y and z
slowuw = [clock' 1*(1+0*sin(2*pi*0.1*clock))'];
slowvw = [clock' 1*(-1+0*sin(2*pi*0.05*clock+pi/5))'];
slowww = [clock' 0+0*sin(2*pi*0.015*clock+pi/3)'];


om_com = 2;

% % Using butter filter
% Wn = om_com*param.dt*2;
% [acc_bcom,acc_acom] = butter(2,Wn);
% Gacc = tf(acc_bcom, acc_acom, param.dt);
% [acc_bcom,acc_acom] = tfdata(d2c(Gacc,'zoh'),'v');
%
% [gyro_bcom,gyro_acom] = butter(2,Wn,'high');
% Ggyro = tf(gyro_bcom, gyro_acom, param.dt);
% [gyro_bcom,gyro_acom] = tfdata(d2c(Ggyro,'zoh'),'v');


% Using simple low-pass filter
acc_bcom = [0 om_com];
acc_acom = [1 om_com];
gyro_bcom = [1 0];
gyro_acom = [1 om_com];
Gacc = tf(acc_bcom,acc_acom);
% [acc_bcom, acc_acom] = tfdata(c2d(Gacc,param.dt),'v');
Ggyro = tf(gyro_bcom,gyro_acom);
% [gyro_bcom, gyro_acom] = tfdata(c2d(Ggyro,param.dt),'v');

bgyro = 1;%[0 10];
agyro = 1;%[1 10];



%
for sigmaw = 0:3 % We design three scenarios: sigma = 0, no wind with noise
    % sigma = 1: no wind psi change
    % sigma = 2: measurement noise with wind
    % sigma = 3: measurement noise with sensor fault, bias angular velocity
    % Twx = [0 10 20 100 110  tfinal];
    % Vwx = [0 0  1   1   0  0]*(sigmaw==1);
    % Twy = [0 10  20  100 110  tfinal];
    % Vwy = [0 0 -0.5 -0.5   0  0]*(sigmaw==1);
    % Vwz = [0 0  0   0   0  0];
    %
    % Vwxsim = [Twx' Vwx'];
    % Vwysim = [Twy' Vwy'];
    % Vwzsim = [Twy' Vwz'];
    
    
    
    
    % for ii = 1:length(clock)
    %     [xref, yref, zref] = quadtrapezoidal(clock(ii));
    %     pm(ii,:) = [xref yref zref];
    % end
    % xrefm = [clock' pm(:,1)];
    % zrefm = [clock' pm(:,3)];
    
    
    % Gain = 0; % If we want to add noise to angular velocity
    %
    % % Simulate noise
    % sigpos = 0.0001;
    % sigeul = 0.00001;
    % sim('quadSimplifiedModel.slx')
    
    
    
    
    % SIMULATION OF ACC AND GYRO NOISES
    sigacc = 1e-6;
    siggyro = 1e-6;
    eacc = sqrt(sigacc)*randn(N,3);
    egyro = sqrt(siggyro)*randn(N,3);
    windConstant = 1;
    windVarying = 1;
    
    % Mass change
    Tm = [0 40 42 120];
    m = [1 1 2 2];
    
    mSim = [Tm' m'];
    
    
    switch(sigmaw)
        case 0 % noise free case
            gainPsi = 0;
            gainWind = 0;
            windConstant = 0;
            windVarying = 0;
            
            Tm = [0 40 42 120];
            m = [1 1 1 1];
            
            mSim = [Tm' m'];
            
            eacc = [clock' eacc*0];
            egyro = [clock' egyro*0];
            sim('quadcopterModel.slx')
            
            stateNon0w0psi = stateNon;
            stateNonSim0w0psi = stateNonsim;
            stateLin0w0psi = stateLin;
            IMU0w0psi = IMU;
            IMUSim0w0psi = IMUSim;
            
            U0w0psi = U;
            pref0w0psi = pref;
            psiref0w0psi = psiref;
            eulref0w0psi = eulref;
            wind0w0psi = wind;
            
            
        case 1 % noise without psi and no wind
            
            gainPsi = 1;
            gainWind = 0;
            windConstant = 0;
            windVarying = 0;
            
            Tm = [0 40 42 120];
            m = [1 1 2 2];
            
            mSim = [Tm' m'];
            
            eacc = [clock' 0*eacc];
            egyro = [clock' 0*egyro];
            sim('quadcopterModel.slx')
            
            stateNon0w1psi = stateNon;
            stateNonSim0w1psi = stateNonsim;
            stateLin0w1psi = stateLin;
            IMU0w1psi = IMU;
            IMUSim0w1psi = IMUSim;
            
            U0w1psi = U;
            pref0w1psi = pref;
            psiref0w1psi = psiref;
            eulref0w1psi = eulref;
            wind0w1psi = wind;
            
        case 2 % contain noise with constant psi and wind
            gainPsi = 0;
            gainWind = 1;
            windConstant = 1;
            windVarying = 1;
            
            eacc = [clock' eacc];
            egyro = [clock' egyro];
            sim('quadcopterModel.slx')
            
            stateNon1w0psi = stateNon;
            stateNonSim1w0psi = stateNonsim;
            stateLin1w0psi = stateLin;
            IMU1w0psi = IMU;
            IMUSim1w0psi = IMUSim;
            
            U1w0psi = U;
            pref1w0psi = pref;
            psiref1w0psi = psiref;
            eulref1w0psi = eulref;
            wind1w0psi = wind;
            
        case 3 % constain noise with change psi and wind
            gainPsi = 1;
            gainWind = 1;
            windConstant = 1;
            windVarying = 1;
            
            eacc = [clock' eacc];
            egyro = [clock' egyro];
            sim('quadcopterModel.slx')
            
            stateNon1w1psi = stateNon;
            stateNonSim1w1psi = stateNonsim;
            stateLin1w1psi = stateLin;
            IMU1w1psi = IMU;
            IMUSim1w1psi = IMUSim;
            
            U1w1psi = U;
            pref1w1psi = pref;
            psiref1w1psi = psiref;
            eulref1w1psi = eulref;
            wind1w1psi = wind;
    end
    % stateNon: pos eul pdot euldot vB pqr vBdot
    % IMU: acc gyro
    % stateLin: p eul pdot euldot peul0w acc wind
    
    % SIMULATION OF LINEAR SYSTEM
    % xref = pref(:,1);
    % yref = pref(:,2);
    % zref = pref(:,3);
    %
    % phiref = eulref(:,1);
    % thetaref = eulref(:,2);
    % psiref = eulref(:,3);
    %
    % x = stateNon(:,1);
    % y = stateNon(:,2);
    % z = stateNon(:,3);
    %
    % phi = stateNon(:,4);
    % theta = stateNon(:,5);
    % psi = stateNon(:,6);
    
    % phi0w = pvector0w(:,4);
    % theta0w = pvector0w(:,5);
    % psi0w = pvector0w(:,6);
    %
    % if sigmaw == 0 % Noise no wind
    %     pos(:,1:3) = [x y z];
    %     gyrom(:,1:3) = gyro;
    %
    %     accm(:,1:3) = acc;
    %     eulrefm(:,1:3) = [phiref thetaref psiref];
    %     eulm(:,1:3) = pvector(:,4:6);
    %     windm(:,1:3) = wind;
    % elseif sigmaw == 1 % Noise and wind
    %     pos(:,4:6) = [x y z];
    %     gyrom(:,4:6) = gyro;
    %     accm(:,4:6) = acc;
    %     eulrefm(:,4:6) = [phiref thetaref psiref];
    %     eulm(:,4:6) = pvector(:,4:6);
    %
    %     windm(:,4:6) = wind;
    % elseif sigmaw == 2 % Noise and no wind and sensor fault
    %     pos0 = [x y z];
    %     gyro0 = gyro;
    %     acc0 = acc;
    %     eulref0 = [phiref thetaref psiref];
    %     eul0 = pvector(:,4:6);
    %     windm(:,7:9) = wind;
    % end
    sigmaw
end

disp('Done simulation')

figPos = [250 300 1100 750];


%% We use Kalman filter to estima phi and theta
% add the bias

% bias_p = [0 kron([0 0 0 1 1 1 0 0 0 -1 -1 -1 0 0 0 1 1 1 0 0 0 0 0 0],ones(1,1000))]';
% bias_q = [0 kron([0 0 0 -1 -1 -1 0 0 0 1 1 1 0 0 0 -1 -1 -1 0 0 0 0 0 0],ones(1,1000))]';
% Gbias = tf(1,[1 2]);
% bias_p = lsim(Gbias,bias_p,clock);
% bias_q = lsim(Gbias,bias_q,clock);


% Ramp up and down
cc0 = 0.0002*(1:2000);
bias_p = [0 zeros(1,2000) cc0 cc0(end)-cc0 -cc0 -cc0(end)+cc0 cc0 cc0(end)-cc0 -cc0 -cc0(end)+cc0,...
    cc0 cc0(end)-cc0 zeros(1,2000)]';
bias_q = -bias_p;

% Ramp
cc0 = [zeros(1,4000) 0.0002*(1:length(clock)-4000)];
bias_p = cc0';
cc0 = [zeros(1,12000)  0.0002*(1:length(clock)-12000)];
bias_q = -bias_p;

% Zero and Constant
cc0 = [zeros(1,4000) 0.0004*(1:1000) 0.4*ones(1,length(clock)-5000)];
bias_p = cc0';
cc0 = [zeros(1,12000) 0.0004*(1:1000) 0.4*ones(1,length(clock)-13000)];
bias_q = -cc0';

% Constant
cc0 = 0.3*ones(N,1);
bias_p = cc0;
bias_q = -cc0;


% Zeros and Constant
cca = 0.5*ones(N,1);
bias_ax = cca;
bias_ay = -cca;

% % Ramp
% cc0 = [zeros(1,6000) 0.0002*(1:length(clock)-6000)];
% bias_ax = cc0';
% cc0 = [zeros(1,10000)  0.0002*(1:length(clock)-10000)];
% bias_ay = -cc0';
% baccLim = 4;

% % trapezoidal
% cc1 = 0.00025*[(1:1000) 1000*ones(1,2000) 1000-(1:1000)];
% bias_ax = [0 zeros(1,2000) cc1 -cc1 cc1 -cc1 zeros(1,6000)]';
% bias_ay = -bias_ax;
% baccLim = 1;

% % sudden bias
zz = zeros(1,2000);
% cc2 = 0.4*ones(1,2000);
% bias_p = [0 zz cc2 zz -cc2 zz cc2 zz -cc2 zz cc2 zz zz]';
% bias_q = -bias_p;
%
% % Sinusoidal
% cc3 = 0.4*sin(2*pi*0.00025*(1:2000));
% bias_ax = [0 zz cc3 -cc3 cc3 -cc3 cc3 -cc3 cc3 -cc3  zz zz zz]';
% bias_ay = -bias_ax;

%%

close all
IMUData = IMU1w0psi;
pnoise = IMUData(:,4);
qnoise = IMUData(:,5);
axnoise = -IMUData(:,1);
aynoise = -IMUData(:,2);
stateData = stateNon1w0psi;
psinoise = stateData(:,6);
psidot = stateData(:,12);

Rpos = 9*0.01; Rspeed = 0.01; Reul = 0.01; Racc  = 1e-6;
xnoise = stateData(:,1) + sqrt(Rpos)*randn(N,1);
ynoise = stateData(:,2) + sqrt(Rpos)*randn(N,1);
unoise = stateData(:,13) + sqrt(Rspeed)*randn(N,1);
vnoise = stateData(:,14) + sqrt(Rspeed)*randn(N,1);
phinoise = stateData(:,4) + sqrt(Reul)*randn(N,1);
thetanoise = stateData(:,5) + sqrt(Reul)*randn(N,1);
index = 200;
% Rpos = 0.01; Rspeed = 0.01; Reul = 0.01; Racc  = 1e-6;

% Tuning coefficient
Qeul = 1e-5*[1 1];
Qv = 1e-5*[1 1];
Qbiasg = 1e-5*[1 1];
Qbiasa = 5*1e-6*[1 1];
Qwind = 1e-4*[1 1];

Qeul = 1e-10*[1 1];
Qv = 1e-10*[1 1];
Qx = 1e-10*[1 1];
Qbiasg = 1e-5*[1 1];
Qbiasa = 1e-6*[1 1];
Qwind = 1e-4*[1 1];



clear x1m y1m x2m y2m x3m x4m x5m x6m y3m y4m y5m y6m
% Define the model, we have two model, one is based on constant mass and
% signal state = [phi theta u v x y biasphi biastheta windx windy]
F1 = [0 0 0 0 0 0 -1 0 0 0;...
    0 0 0 0 0 0 0 -1 0 0;...
    0 param.g(3), -param.lambda1/param.m 0, 0 0,0 0, -param.lambda1/param.m 0;...
    -param.g(3) 0, 0 -param.lambda1/param.m, 0 0,0 0, 0 -param.lambda1/param.m;...
    0 0 1 0 0 0 0 0 0 0;...
    0 0 0 1 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0];F1 = eye(10) + F1*param.dt;
B1 = [param.dt 0 0 0 0 0 0 0 0 0;0 param.dt 0 0 0 0 0 0 0 0]';
G1 = [1 0 0 0 0 0 0 0;...
    0 1 0 0 0 0 0 0;...
    0 0 param.dt 0 0 0 0 0;...
    0 0 0 param.dt 0 0 0 0;...
    0 0 param.dt^2/2 0 0 0 0 0;...
    0 0 0 param.dt^2/2 0 0 0 0;...
    0 0 0 0 1 0 0 0;...
    0 0 0 0 0 1 0 0;...
    0 0 0 0 0 0 1 0;...
    0 0 0 0 0 0 0 1];
H11 = [0 0 0 0 1 0 0 0 0 0;...
    0 0 0 0 0 1 0 0 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0];
H12 = [0 0 param.lambda1/param.m 0 0 0 0 0 param.lambda1/param.m 0;...
    0 0 0 param.lambda1/param.m 0 0 0 0 0 param.lambda1/param.m];

H1 = [H11;H12];
R11 = diag([Rpos Rpos Rspeed Rspeed]);
R12 = diag([Racc Racc ]);
R1 = blkdiag(R11,R12);
Q1 = diag([Qeul, Qv, Qbiasg, Qwind]);


u1 = [pnoise+bias_p qnoise+bias_q];
y1 = [xnoise ynoise unoise vnoise axnoise aynoise];

% Stationary Kalman
% [K1bar,Pp1,Pf1] = dlqe(F1,G1,[H11;H12],Q1,R1);
% [K11bar,Pp11,Pf11] = dlqe(F1,G1,H11,Q1,R11);

x1 = [zeros(8,1);1;-1];
x1m =  zeros(N,10);
y1m = zeros(N,6);
P1 = 1000*eye(10);
for ii = 1:length(clock)
    % Prediction step
    x1 = F1*x1+B1*u1(ii,:)';
    P1 = F1*P1*F1'+G1*Q1*G1';
    
    % Correction step
    if mod(ii-1,index) == 0
        S1 = H1*P1*H1'+R1;
        K1 = P1*H1'/(S1);
        x1 = x1 + K1*(y1(ii,:)' - H1*x1);
        P1 = P1 - K1*H1*P1;
    else
        S1 = H12*P1*H12'+R12;
        K1 = P1*H12'/(S1);
        x1 = x1 + K1*(y1(ii,5:end)' - H12*x1);
        P1 = P1 - K1*H12*P1;
    end
    
    x1m(ii,:) = x1;
    y1m(ii,:) = H1*x1;
    
end

% The second approach, taking ay as the input as well
% x = [phi theta u v x y biasphi biastheta]
F2 = [0 0 0 0 0 0 -1 0;...
    0 0 0 0 0 0 0 -1;...
    0 param.g(3) 0 0 0 0 0 0;...
    -param.g(3) 0 0 0 0 0 0 0;...
    0 0 1 0 0 0 0 0;...
    0 0 0 1 0 0 0 0;...
    0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0]; F2 = eye(8)+F2*param.dt;
B2 = param.dt*[1 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0; 0 0 -1 0 0 0 0 0;0 0 0 -1 0 0 0 0]';
G2 = [1 0 0 0 0 0;...
    0 1 0 0 0 0;...
    0 0 param.dt 0 0 0;...
    0 0 0 param.dt 0 0;...
    0 0 param.dt^2/2 0 0 0;...
    0 0 0 param.dt^2/2 0 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 1];
H21 = [0 0 0 0 1 0 0 0;...
    0 0 0 0 0 1 0 0;...
    0 0 1 0 0 0 0 0;...
    0 0 0 1 0 0 0 0];
H22 = [];
H2 = [H21;H22];

u2 = [pnoise+bias_p qnoise+bias_q axnoise aynoise];
y2 = [xnoise ynoise unoise vnoise];

R21 = diag([Rpos Rpos Rspeed Rspeed]);
R22 = diag([]);
R2 = blkdiag(R21,R22);
Q2 = diag([Qeul, Qv, Qbiasg]);

% Stationary Kalman
[K2bar,Pp2,Pf2] = dlqe(F2,G2,H2,Q2,R2);

x2 = zeros(8,1);
x2m = zeros(N,8);
y2m = zeros(N,4);
P2 = 1000*eye(8);
for ii = 1:length(clock)
    % Prediction step
    x2 = F2*x2+B2*u2(ii,:)';
    P2 = F2*P2*F2'+G2*Q2*G2';
    
    if mod(ii-1,index) == 0
        % Correction step
        S2 = H2*P2*H2'+R2;
        K2 = P2*H2'/(S2);
        x2 = x2+K2*(y2(ii,:)'-H2*x2);
        P2 = P2 - K2*H2*P2;
    else
        %        % Correction step
        %        S2 = H22*P2*H22'+R22;
        %        K2 = P2*H22'*inv(S2);
        %        x2 = x2+K2*(y2(ii,3:end)'-H22*x2);
        %        P2 = P2 - K2*H22*P2;
    end
    
    
    x2m(ii,:) = x2;
    y2m(ii,:) = H2*x2;
end


% The second approach where bias comes from the acc measurement
% x = [phi theta u v x y biasx biasy windx windy]
F3 = [  0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 param.g(3), -param.lambda1/param.m 0, 0 0, 0 0, -param.lambda1/param.m 0;...
    -param.g(3) 0, 0 -param.lambda1/param.m, 0 0, 0 0, 0 -param.lambda1/param.m;...
    0 0 1 0 0 0 0 0 0 0;...
    0 0 0 1 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0];F3 = eye(10)+param.dt*F3;

B3 = param.dt*[1 0 0 0 0 0 0 0 0 0;...
    0 1 0 0 0 0 0 0 0 0]';
H31 = [0 0 0 0 1 0 0 0 0 0;...
    0 0 0 0 0 1 0 0 0 0;...
    0 0 1 0 0 0 0 0 0 0;...
    0 0 0 1 0 0 0 0 0 0];
H32 = [0 0 param.lambda1/param.m 0 0 0 1 0 param.lambda1/param.m 0;...
    0 0 0 param.lambda1/param.m 0 0 0 1 0 param.lambda1/param.m];
H3 = [H31;H32];
G3 = [1 0 0 0 0 0 0 0;...
    0 1 0 0 0 0 0 0;...
    0 0 param.dt 0 0 0 0 0;...
    0 0 0 param.dt 0 0 0 0;...
    0 0 param.dt^2/2 0 0 0 0 0;...
    0 0 0 param.dt^2/2 0 0 0 0;...
    0 0 0 0 1 0 0 0;...
    0 0 0 0 0 1 0 0;...
    0 0 0 0 0 0 1 0;...
    0 0 0 0 0 0 0 1];

u3 = [pnoise qnoise];
y3 = [xnoise ynoise unoise vnoise axnoise+bias_ax aynoise+bias_ay];

R31 = diag([Rpos Rpos Rspeed Rspeed]);
R32 = diag([Racc Racc]);
R3 = blkdiag(R31,R32);
Q3 = diag([Qeul, Qv, Qbiasa, Qwind]);
% Stationary Kalman
[K3bar,Pp3,Pf3] = dlqe(F3,G3,H3,Q3,R3);

x3 = [zeros(8,1);1;-1];
x3m = zeros(N,10);
y3m = zeros(N,6);
P3 = 1000*eye(10);
for ii = 1:length(clock)
    % Prediction step
    x3 = F3*x3 + B3*u3(ii,:)';
    P3 = F3*P3*F3' + G3*Q3*G3';
    
    if mod(ii-1,index) == 0
        % Correction step
        S3 = H3*P3*H3'+blkdiag(R31,R32);
        K3 = P3*H3'/(S3);
        x3 = x3 + K3*(y3(ii,:)' - H3*x3);
        P3 = P3 - K3*H3*P3;
    else
        % Correction step
        S3 = H32*P3*H32'+R32;
        K3 = P3*H32'/(S3);
        x3 = x3 + K3*(y3(ii,5:end)' - H32*x3);
        P3 = P3 - K3*H32*P3;
    end
    
    x3m(ii,:) = x3;
    y3m(ii,:) = H3*x3;
end



% Using wind estimate as well as bias estimate
% x = [phi theta u v x y biasax biasay]
F4 = [0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0;...
    0 param.g(3) 0 0 0 0 1 0;...
    -param.g(3) 0 0 0 0 0 0 1;...
    0 0 1 0 0 0 0 0;...
    0 0 0 1 0 0 0 0;...
    0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0]; F4 = eye(8)+F4*param.dt;
H41 = [0 0 0 0 1 0 0 0;...
    0 0 0 0 0 1 0 0;...
    0 0 1 0 0 0 0 0;...
    0 0 0 1 0 0 0 0];

H42 = [];
H4 = [H41;H42];
B4 = param.dt*[1 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0; 0 0 -1 0 0 0 0 0;0 0 0 -1 0 0 0 0]';
G4 = [1 0 0 0 0 0;...
    0 1 0 0 0 0;...
    0 0 param.dt 0 0 0;...
    0 0 0 param.dt 0 0;...
    0 0 param.dt^2/2 0 0 0;...
    0 0 0 param.dt^2/2 0 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 1];

u4 = [pnoise qnoise axnoise+bias_ax aynoise+bias_ay ];
y4 = [xnoise ynoise unoise vnoise];

R41 = diag([Rpos Rpos Rspeed Rspeed]);
R42 =  diag([]);
R4 = blkdiag(R41,R42);
Q4 = diag([Qeul, Qv, Qbiasa]);
% Stationary Kalman
[K4bar,Pp4,Pf4] = dlqe(F4,G4,H4,Q4,R4);

x4 = zeros(8,1);
x4m = zeros(N,8);
y4m = zeros(N,4);
P4 = 1000*eye(8);
for ii = 1:length(clock)
    % Prediction step
    x4 = F4*x4 + B4*u4(ii,:)';
    P4 = F4*P4*F4' + G4*Q4*G4';
    
    if mod(ii-1,index) == 0
        % Correction step
        S4 = H4*P4*H4'+R4;
        K4 = P4*H4'/(S4);
        x4 = x4 + K4*(y4(ii,:)' - H4*x4);
        P4 = P4 - K4*H4*P4;
        
    else
        %        % Correction step
        %        S4 = H42*P4*H42'+R42;
        %        K4 = P4*H42'*inv(S4);
        %        x4 = x4 + K4*(y4(ii,3:end)' - H42*x4);
        %        P4 = P4 - K4*H42*P4;
    end
    x4m(ii,:) = x4;
    y4m(ii,:) = H4*x4;
end


% ESTIMATE BOTH GYRO BIAS AND ACC BIAS
% Using wind estimate as well as bias estimate
% x = [phi theta u v x y biasphi biastheta biasx biasy windx windy]
F5 = [  0 0, 0 0, 0 0, -1 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 -1, 0 0, 0 0;...
    0 param.g(3), -param.lambda1/param.m 0, 0 0, 0 0, 0 0, -param.lambda1/param.m 0;...
    -param.g(3) 0, 0 -param.lambda1/param.m, 0 0, 0 0, 0 0, 0 -param.lambda1/param.m;...
    0 0, 1 0, 0 0, 0 0, 0 0, 0 0;...
    0 0 0 1 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0];F5 = eye(12)+param.dt*F5;
H51 = [0 0, 0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0, 0 0];
H52 = [0 0, param.lambda1/param.m 0, 0 0, 0 0, 1 0, param.lambda1/param.m 0;...
    0 0, 0 param.lambda1/param.m, 0 0, 0 0, 0 1, 0 param.lambda1/param.m];
H5 = [H51;H52];
B5 = param.dt*[1 0 0 0 0 0 0 0 0 0 0 0;...
    0 1 0 0 0 0 0 0 0 0 0 0]';
G5 = [1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0, 0 0;...
    0 0, param.dt^2/2 0, 0 0, 0 0, 0 0;...
    0 0, 0 param.dt^2/2, 0 0, 0 0, 0 0;...
    0 0, 0 0, 1 0 0 0 0 0;...
    0 0, 0 0, 0 1 0 0 0 0;...
    0 0, 0 0, 0 0 1 0 0 0;...
    0 0, 0 0, 0 0 0 1 0 0;...
    0 0, 0 0, 0 0 0 0, 1 0;...
    0 0, 0 0, 0 0 0 0, 0 1];

u5 = [pnoise+bias_p qnoise+bias_q];
y5 = [xnoise ynoise unoise vnoise axnoise+bias_ax aynoise+bias_ay];

R51 = diag([Rpos Rpos Rspeed Rspeed]);
R52 = diag([Racc Racc]);
R5 = blkdiag(R51,R52);
Q5 = diag([Qeul, Qv, Qbiasg, Qbiasa, Qwind]);
% Stationary Kalman
[K5bar,Pp5,Pf5] = dlqe(F5,G5,H5,Q5,R5);

x5 = [zeros(10,1);1 ;-1];
x5m = zeros(N,12);
y5m = zeros(N,6);
P5 = 1000*eye(12);
for ii = 1:length(clock)
    % Prediction step
    x5 = F5*x5 + B5*u5(ii,:)';
    P5 = F5*P5*F5' + G5*Q5*G5';
    
    if mod(ii-1,index) == 0
        % Correction step
        S5 = H5*P5*H5'+R5;
        K5 = P5*H5'/(S5);
        x5 = x5 + K5*(y5(ii,:)' - H5*x5);
        P5 = P5 - K5*H5*P5;
    else
        % Correction step
        S5 = H52*P5*H52'+R52;
        K5 = P5*H52'/(S5);
        x5 = x5 + K5*(y5(ii,5:end)' - H52*x5);
        P5 = P5 - K5*H52*P5;
    end
    
    x5m(ii,:) = x5;
    y5m(ii,:) = H5*x5;
end

% USING AX AY AS THE INPUT
% x = [phi theta u v x y biasphi biastheta biasx biasy]
F6 = [  0 0, 0 0, 0 0, -1 0, 0 0;...
    0 0, 0 0, 0 0, 0 -1, 0 0;...
    0 param.g(3), 0 0, 0 0, 0 0, 1 0;...
    -param.g(3) 0, 0 0, 0 0, 0 0, 0 1;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];F6 = eye(10)+param.dt*F6;
H61 = [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0];
H62 = [];
H6 = [H61;H62];
B6 = param.dt*[1 0 0 0 0 0 0 0 0 0;...
    0 1 0 0 0 0 0 0 0 0;...
    0 0, -1 0, 0 0, 0 0,0 0;...
    0 0, 0 -1, 0 0, 0 0, 0 0]';
G6 = [1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0;...
    0 0, param.dt^2/2 0, 0 0, 0 0;...
    0 0, 0 param.dt^2/2, 0 0, 0 0;...
    0 0, 0 0, 1 0 0 0;...
    0 0, 0 0, 0 1 0 0;...
    0 0, 0 0, 0 0 1 0;...
    0 0, 0 0, 0 0 0 1];

u6 = [pnoise+bias_p qnoise+bias_q axnoise+bias_ax aynoise+bias_ay];
y6 = [xnoise ynoise unoise vnoise];

R61 = diag([Rpos Rpos Rspeed Rspeed]);...
    R62 = diag([]);
R6 = blkdiag(R61,R62);
Q6 = 10*diag([Qeul, Qv, Qbiasg, Qbiasa]);
% Stationary Kalman
[K6bar,Pp6,Pf6] = dlqe(F6,G6,H6,Q6,R6);

x6 = zeros(10,1);
x6m = zeros(N,10);
y6m = zeros(N,4);
P6 = 1000*eye(10);
for ii = 1:length(clock)
    % Prediction step
    x6 = F6*x6 + B6*u6(ii,:)';
    P6 = F6*P6*F6' + G6*Q6*G6';
    
    if mod(ii-1,index) == 0
        % Correction step
        S6 = H6*P6*H6'+R6;
        K6 = P6*H6'/S6;
        x6 = x6 + K6*(y6(ii,:)' - H6*x6);
        P6 = P6 - K6*H6*P6;
    else
        %        % Correction step
        %        S6 = H62*P6*H62'+R62;
        %        K6 = P6*H62'*inv(S6);
        %        x6 = x6 + K6*(y6(ii,3:end)' - H62*x6);
        %        P6 = P6 - K6*H62*P6;
    end
    
    x6m(ii,:) = x6;
    y6m(ii,:) = H6*x6;
end
disp('Done consider Constant Yaw angle')

%

x5Smooth = x5m;
x6Smooth = x6m;
for ii = 1:N
    if mod(ii-1,index) == 0
        if ii<N-index+1
            rate = (x6m(ii+index,9:10)-x6m(ii,9:10))/index;
            x6Smooth(ii+1:ii+index-1,9:10) = x6Smooth(ii,9:10)+rate.*(1:index-1)';
            rate = (x5m(ii+index,9:10)-x5m(ii,9:10))/index;
            x5Smooth(ii+1:ii+index-1,9:10) = x5Smooth(ii,9:10)+rate.*(1:index-1)';
        end
    end
end
%

beulLim = 0.6;horix = -1900 ;very = 700;movex = 450; movey = 450; figx = 420; figy = 320;%560 420
% -1950 700 420 320 450 300
close all
ff(1) = figure(1);
ff(1).Position = [horix very figx figy];
plot(clock,stateNon1w0psi(:,7),'linewidth',1.5)
hold on
plot(clock,x1m(:,3),'--','linewidth',1.5)
plot(clock,x2m(:,3),'--','linewidth',1.5)
legend 'u true' 'u GYRO' 'u GYRO accin'% 'u ax bias'

ff(2) = figure(2);
ff(2).Position = [horix+movex very figx figy];
h2 = axes;h2.FontSize = 14;hold on
subplot(2,1,1);
plot(clock,bias_p,'linewidth',1.5)
hold on
plot(clock,x1m(:,7),'--','linewidth',1.5)
plot(clock,x2m(:,7),'-.','linewidth',1.5)
legend 'b-phi true' 'b-phi Aout GYRO' 'b-phi Ain GYRO'% 'v ax bias'
xlabel 'Time [s]'
ylim([-beulLim beulLim])
grid on

subplot(2,1,2)
plot(clock,bias_q,'linewidth',1.5)
hold on
plot(clock,x1m(:,8),'--','linewidth',1.5)
plot(clock,x2m(:,8),'-.','linewidth',1.5)
legend 'b-theta true' 'b-theta Aout GYRO' 'b-theta Ain GYRO'% 'v ax bias'
xlabel 'Time [s]'
ylim([-beulLim beulLim])
grid on


ff(3) = figure(3);
ff(3).Position = [horix+2*movex very figx figy];
h3 = axes;h3.FontSize = 14;hold on
plot(clock,wind1w0psi(:,1),'linewidth',1.5)
hold on
plot(clock,x1m(:,9),'linewidth',1.5)
plot(clock,x3m(:,9),'linewidth',1.5)
plot(clock,x5m(:,11),'linewidth',1.5)
legend 'u_w true' 'u_w GYRO' 'u_w Acc' 'u_w Gyro+acc'
xlabel 'Time [s]'
ylabel 'u_w [m/s]'
ylim([-4 4])

ff(4) = figure(4);
ff(4).Position = [horix+3*movex very figx figy];
h4 = axes;h4.FontSize = 14;hold on
plot(clock,wind1w0psi(:,2),'linewidth',1.5)
hold on
plot(clock,x1m(:,10),'linewidth',1.5)
plot(clock,x3m(:,10),'linewidth',1.5)
plot(clock,x5m(:,12),'linewidth',1.5)
legend 'v_w true' 'v_w GYRO' 'v_w Acc' 'v_w Gyro+acc'
xlabel 'Time [s]'
ylabel 'v_w [m/s]'
ylim([-4 4])

%
ff(5) = figure(5);clf % bias ax
ff(5).Position = [horix very-movey figx figy];
subplot(2,1,1)
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,x3m(:,7),'--','linewidth',1.5)
plot(clock,x4m(:,7),'-.','linewidth',1.5)
ylim([-baccLim baccLim])
legend 'b-a_x' 'b-a_x Aout Acc' 'b-a_x Ain Acc'
xlabel 'Time [s]'
grid on


subplot(2,1,2)
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,x3m(:,8),'--','linewidth',1.5)
plot(clock,x4m(:,8),'-.','linewidth',1.5)
ylim([-baccLim baccLim])
legend 'b-a_y' 'b-a_y Aout Acc' 'b-a_y Aout Acc'
xlabel 'Time [s]'
grid on


ff(6) = figure(6);clf % bias ay
ff(6).Position = [horix+movex very-movey figx figy];
subplot(2,1,1)
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,x5m(:,9),'-','linewidth',1.5)
plot(clock,x6Smooth(:,9),'-.','linewidth',1.5)
ylim([-baccLim baccLim])
legend 'b-a_x' 'b-a_x Aout Gyro+Acc' 'b-a_x Ain Gyro+Acc'
xlabel 'Time [s]'
grid on


subplot(2,1,2)
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,x5m(:,10),'-','linewidth',1.5)
plot(clock,x6Smooth(:,10),'-.','linewidth',1.5)
ylim([-baccLim baccLim])
legend 'b-a_y' 'b-a_y Aout Gyro+Acc' 'b-a_y Aout Gyro+Acc'
xlabel 'Time [s]'
grid on


%
factor = 1;
ff(7) = figure(7);clf % bias p
ff(7).Position = [horix+2*movex very-movey figx figy];
h7 = axes;h7.FontSize = 14;hold on
subplot(2,1,1)
plot(clock(1:factor:end),bias_p(1:factor:end),'linewidth',1.5)
hold on
plot(clock(1:factor:end),x5m(1:factor:end,7),'--','linewidth',1.5)
plot(clock(1:factor:end),x6m(1:factor:end,7),'-.','linewidth',1.5)
ylim([-beulLim beulLim])
xlabel 'Time [s]'
ylabel '\beta_\phi'
legend 'b-phi true' 'b-phi Aout GYRO+Acc' 'b-phi Ain GYRO+Acc'
grid on

subplot(2,1,2)
plot(clock(1:factor:end),bias_q(1:factor:end),'linewidth',1.5)
hold on
plot(clock(1:factor:end),x5m(1:factor:end,8),'--','linewidth',1.5)
plot(clock(1:factor:end),x6m(1:factor:end,8),'-.','linewidth',1.5)
ylim([-beulLim beulLim])
xlabel 'Time [s]'
ylabel '\beta_\theta'
legend 'b-theta true' 'b-theta Aout GYRO+Acc' 'b-theta Ain GYRO+Acc'
grid on


ff(8) = figure(8); % bias q
ff(8).Position = [horix+3*movex very-movey figx figy];
h8 = axes;h8.FontSize = 14;hold on
hold on

[rank(obsv(F1,H1)),rank(obsv(F2,H2)),rank(obsv(F3,H3)),rank(obsv(F4,H4)),...
    rank(obsv(F5,H5)),rank(obsv(F6,H6))]

%%
% Now we consider the sensor-to-sensor model for mass change detection
% param.m = 0.5;
Gs = -tf([0 param.g(3)*param.lambda1/param.m],[1 param.lambda1/param.m 0]);
Gd = c2d(Gs,dt,'tustin');
[bhat,ahat] = tfdata(Gd,'v');
th = [ahat(2:end) (bhat(1:end))]';

Nw = 500;
lambdaW = ones(Nw,1);
for ii = 1:Nw
    lambdaW(ii,:) = 0.999^(Nw-ii);
end
lambdaW = sqrt(lambdaW);

% no psi case
pnoise = -IMU1w0psi(:,4);
qnoise = IMU1w0psi(:,5);
axnoise = IMU1w0psi(:,1);
aynoise = IMU1w0psi(:,2);

q0 = IMU0w0psi(:,5);
ax0 = IMU0w0psi(:,1);
p0 = -IMU0w0psi(:,4);
ay0 = IMU0w0psi(:,2);

% Psi case
phatpsi = -IMU1w1psi(:,4);
qhatpsi = IMU1w1psi(:,5);
axpsi = IMU1w1psi(:,1);
aypsi = IMU1w1psi(:,2);

q0psi = IMU0w1psi(:,5);
ax0psi = IMU0w1psi(:,1);
p0psi = -IMU0w1psi(:,4);
ay0psi = IMU0w1psi(:,2);

qhat = qnoise + (bias_q - x6m(:,8));
phat =  pnoise + (bias_p - x6m(:,7));
ay = aynoise + (bias_ay - x6Smooth(:,10));
ax = axnoise + (bias_ax - x6Smooth(:,9));

Zx0 = [-ax0(2:end-1) -ax0(1:end-2) q0(3:end) q0(2:end-1) q0(1:end-2)];
Phix = [-ax(2:end-1) -ax(1:end-2) qhat(3:end) qhat(2:end-1) qhat(1:end-2)];
Yx = ax(3:end);

Zy0 = [-ay0(2:end-1) -ay0(1:end-2) p0(3:end) p0(2:end-1) p0(1:end-2)];
Phiy = [-ay(2:end-1) -ay(1:end-2) phat(3:end) phat(2:end-1) phat(1:end-2)];
Yy = ay(3:end);

Zx0psi = [-ax0psi(2:end-1) -ax0psi(1:end-2) q0psi(3:end) q0psi(2:end-1) q0psi(1:end-2)];
Phixpsi = [-axpsi(2:end-1) -axpsi(1:end-2) qhatpsi(3:end) qhatpsi(2:end-1) qhatpsi(1:end-2)];
Yxpsi = axpsi(3:end);

Zy0psi = [-ay0(2:end-1) -ay0(1:end-2) p0(3:end) p0(2:end-1) p0(1:end-2)];
Phiypsi = [-aypsi(2:end-1) -aypsi(1:end-2) phatpsi(3:end) phatpsi(2:end-1) phatpsi(1:end-2)];
Yypsi = aypsi(3:end);

Nn = size(Zy0,1);
for ii = Nw:Nn
    Z0yii = [Zy0(ii-Nw+1:ii,:) [zeros(1,5);Zy0(ii-Nw+1:ii-1,:)] [zeros(2,5);Zy0(ii-Nw+1:ii-2,:)],...
        [zeros(3,5);Zy0(ii-Nw+1:ii-3,:)] [zeros(4,5);Zy0(ii-Nw+1:ii-4,:)] [zeros(5,5);Zy0(ii-Nw+1:ii-5,:)]];
    Phiyii = Phiy(ii-Nw+1:ii,:);
    Yyii = Yy(ii-Nw+1:ii,:);
    Eyii = Yyii-Phiyii*th;
    YY = norm(Z0yii'*(Eyii.*lambdaW));
    Jyii(ii,:) = YY;
    Jyii2(ii,:) = norm(Eyii);
    
    Z0xii = Zx0(ii-Nw+1:ii,:);
    Phixii = Phix(ii-Nw+1:ii,:);
    Yxii = Yx(ii-Nw+1:ii,:);
    Exii = Yxii-Phixii*th;
    XX = norm(Z0xii'*(Exii.*lambdaW));
    Jxii(ii,:) = XX;
    
    Z0yiipsi = [Zy0psi(ii-Nw+1:ii,:)];
    Phiyiipsi = Phiypsi(ii-Nw+1:ii,:);
    Yyiipsi = Yypsi(ii-Nw+1:ii,:);
    Eyiipsi = Yyiipsi-Phiyiipsi*th;
    YYpsi = norm(Z0yiipsi'*(Eyiipsi.*lambdaW));
    Jyiipsi(ii,:) = YYpsi;
    
    Z0xiipsi = Zx0psi(ii-Nw+1:ii,:);
    Phixiipsi = Phixpsi(ii-Nw+1:ii,:);
    Yxiipsi = Yxpsi(ii-Nw+1:ii,:);
    Exiipsi = Yxiipsi-Phixiipsi*th;
    XXpsi = norm(Z0xiipsi'*(Exiipsi.*lambdaW));
    Jxiipsi(ii,:) = XXpsi;
end
%
% IF extended IV, normalize
basicIV = mean(Jxii(12000:end));
extendIV = mean(Jyii(12000:end));


%
h9 = figure;clf
h9.Position =  [horix very figx figy];
h99 = axes;h99.FontSize = 14;hold on
plot(clock(1,1:Nn),Jyii/extendIV*basicIV)
hold on
plot(clock(1,1:Nn),Jxii)
% plot(clock(1,1:Nn),Jyii2)
% plot(clock,bias_p,'linewidth',2)
% plot(clock,bias_q,'linewidth',2)
% plot(clock(1,1:Nn),Jyiipsi)
% plot(clock(1,1:Nn),Jxiipsi)
xlabel 'Time [s]'
ylabel 'IV cost function J_{IV}'
legend 'Pitch model, Extended IV' 'Roll model, basic IV'% 'pitch psi' 'roll psi'
grid on
ylim([-0.1 0.3])
% h99.YTick = [0 0.05 0.1 0.15];
% ylim([0 1])
















%%
% stateNon: pos eul pdot euldot vB pqr vBdot
% IMU: acc gyro
% stateLin: p eul pdot euldot peul0w acc wind
% GPS sampling time, 200 means 200x0.005 = 1 second
index = 200;
% Taking the measurements
Rpos = 9*0.01; Rspeed = 0.01; Reul = 0.01; Racc  = 1e-6;

IMUdata = IMU1w1psi;
stateData = stateNon1w1psi;
% IMUData = [IMU0w1psi(:,1:3) stateNon0w1psi(:,19:21)];

windData = wind1w1psi;

pnoise = IMUdata(:,4);
qnoise = IMUdata(:,5);
rnoise = IMUdata(:,6);
axnoise = -IMUdata(:,1);
aynoise = -IMUdata(:,2);

xnoise = stateData(:,1) + 0*sqrt(Rpos)*randn(N,1);
ynoise = stateData(:,2) + 0*sqrt(Rpos)*randn(N,1);
xdotnoise = stateData(:,7) + 0*sqrt(Rspeed)*randn(N,1);
ydotnoise = stateData(:,8) + 0*sqrt(Rspeed)*randn(N,1);

unoise = stateData(:,13) + 0*sqrt(Rspeed)*randn(N,1);
vnoise = stateData(:,14) + 0*sqrt(Rspeed)*randn(N,1);
wnoise = stateData(:,15);
phinoise = stateData(:,4) + 0*sqrt(Reul)*randn(N,1);
thetanoise = stateData(:,5) + 0*sqrt(Reul)*randn(N,1);
udot = stateData(:,19);
vdot = stateData(:,20);

psinoise = stateData(1:end,6);
psidot = stateData(1:end,18); % Similar performance to state 18, r
psidot2 = stateData(:,18) + stateData(:,4).*stateData(:,4);
% Check the model, x = [phi theta u v x y]
p = [0 0]';p2 = [0 0]';p3 = [0 0]';p0 = p;
eulSim = [0;0];vecSim = [0;0];vecCal = [0;0];
for ii = 1:length(xnoise)
    % Sim Euler angle
    Feul = [0 psidot(ii);-psidot(ii) 0];
    eulSim = (eye(2) + param.dt*Feul)*eulSim + param.dt*[pnoise(ii);qnoise(ii)];
    
    % Sim velocity and position
    Fvec = [0 param.g(3) -param.lambda1/param.m psidot(ii);...
        -param.g(3) 0 -psidot(ii) -param.lambda1/param.m];
    Fpos = [cos(psinoise(ii)) -sin(psinoise(ii));...
        sin(psinoise(ii)) cos(psinoise(ii))];
    
    vecdot = Fvec*[phinoise(ii);thetanoise(ii);vecSim]-param.lambda1/param.m*Fpos'*windData(ii,1:2)';
    vecSim = vecSim+param.dt*vecdot;
    
    posdot = Fpos*vecSim;
    posdottrue = Fpos*[unoise(ii);vnoise(ii)];
    
    IMUatrue = param.lambda1/param.m*eye(2)*[unoise(ii);vnoise(ii)] + param.lambda1/param.m*Fpos'*windData(ii,1:2)';
    vecdotCal = param.g(3)*[eulSim(2);-eulSim(1)]+Feul*[unoise(ii);vnoise(ii)]-[axnoise(ii);aynoise(ii)];
    vecCal = vecCal + param.dt*vecdotCal;
    
    Rii = quadeul2rot([phinoise(ii);thetanoise(ii);psinoise(ii)]);
    pdot = Fpos*[unoise(ii);vnoise(ii)];
    pdot2 = stateData(ii,7:8)';
    pdot3 = Rii*stateData(ii,13:15)';
    
    IMUatruem(ii,:) = IMUatrue;
    vecdotCalm(ii,:) = vecdotCal;
    
    eulSimm(ii,:) = eulSim;
    
    vecSimm(ii,:) = vecSim;
    vecCalm(ii,:) = vecCal;
    vecdotm(ii,:) = vecdot;
    
    posdottruem(ii,:) = posdottrue;
    p0m(ii,:) = p0;
    pm(ii,:) = p;
    pm2(ii,:) = p2;
    pm3(ii,:) = p3;
    p0 = p0 + param.dt*posdot;
    p = p + param.dt*pdot;
    p2 = p2 + param.dt*pdot2;
    p3 = p3 + param.dt*pdot3(1:2);
end
unoise = stateData(:,13) + 0*sqrt(Rspeed)*randn(N,1);
vnoise = stateData(:,14) + 0*sqrt(Rspeed)*randn(N,1);
xdotnoise = stateData(:,7) + 0*sqrt(Rspeed)*randn(N,1);
ydotnoise = stateData(:,8) + 0*sqrt(Rspeed)*randn(N,1);
%
close all
fff(1) = figure(1);clf
fff(1).Position = [-1950 700 560 420];
subplot(2,1,1)
plot(eulSimm(:,1),'linewidth',1.5)
hold on
plot(phinoise,'--','linewidth',1.5)
legend 'phi cal' 'phi true'
subplot(2,1,2)
plot(eulSimm(:,2),'linewidth',1.5)
hold on
plot(thetanoise,'--','linewidth',1.5)
legend 'theta cal' 'theta true'

fff(2) = figure(2);clf
fff(2).Position = [-1950+600 700 560 420];
plot(xnoise,'linewidth',1.5)
hold on
plot(p0m(:,1),'--','linewidth',1.5)
plot(pm(:,1),'--','linewidth',1.5)
plot(pm2(:,1),'-.','linewidth',1.5)
plot(pm3(:,1),'--','linewidth',1.5)
legend 'true x' 'cal from int a' 'cal from uv*R' 'cal from xydot' 'cal from xyzdot*R'

fff(3) = figure(3);clf
fff(3).Position = [-1950+2*600 700 560 420];
subplot(2,1,1)
plot(unoise,'linewidth',1.5)
hold on
plot(vecSimm(:,1),'-.','linewidth',1.5)
plot(vecCalm(:,1),'-.','linewidth',1.5)
legend 'u true' 'u int a' 'u from a'


subplot(2,1,2)
plot(vnoise,'linewidth',1.5)
hold on
plot(vecSimm(:,2),'-.','linewidth',1.5)
plot(vecCalm(:,2),'-.','linewidth',1.5)
legend 'v true' 'v int a' 'v from a'

fff(4) = figure(4);clf
fff(4).Position = [-1950 300 560 420];
subplot(2,1,1)
plot(xdotnoise,'linewidth',1.5)
hold on
plot(posdottruem(:,1),'-.','linewidth',1.5)


subplot(2,1,2)
plot(ydotnoise,'linewidth',1.5)
hold on
plot(posdottruem(:,2),'-.','linewidth',1.5)

fff(5) = figure(5);clf
fff(5).Position = [-1950+600 300 560 420];
subplot(2,1,1)
plot(udot,'linewidth',1.5)
hold on
plot(vecdotCalm(:,1),'-.','linewidth',1.5)
legend 'udot'

subplot(2,1,2)
plot(vdot,'linewidth',1.5)
hold on
plot(vecdotCalm(:,2),'-.','linewidth',1.5)

legend 'vdot'

fff(6) = figure(6);clf
fff(6).Position = [-1950+2*600 300 560 420];
subplot(2,1,1)
plot(axnoise,'linewidth',1.5)
hold on
plot(IMUatruem(:,1),'-.','linewidth',1.5)
legend 'udot'

subplot(2,1,2)
plot(aynoise,'linewidth',1.5)
hold on
plot(IMUatruem(:,2),'-.','linewidth',1.5)

legend 'vdot'


%%
% Tuning coefficient, nonlinear function should be consider where the rate
% goes to the system in the euler angle as well as acceleration
Qeul = 1e-10*[1 1];
Qv = 1e-10*[1 1];
Qx = 1e-10*[1 1];
Qbiasg = 1e-5*[1 1];
Qbiasa = 1e-6*[1 1];
Qwind = 1e-4*[1 1];

% Qeul = 1e-10*[1 1];
% Qv = 1e-5*[1 1];
% Qx = 1e-5*[1 1];
% Qbiasg = 1e-5*[1 1];
% Qbiasa = 1e-6*[1 1];
% Qwind = 1e-5*[1 1];

% Estimate gyro bias without euler angle measurement, we use yaw
% angle change instead of
% x = [phi theta u v x y biasphi biastheta windx windy]
clear x5m y5m x6m y6m P5m P6m
F5 = @(psit, r) eye(10)+param.dt*[  0 r, 0 0, 0 0, -1 0, 0 0;...
    -r 0, 0 0, 0 0, 0 -1, 0 0;...
    0 param.g(3), -param.lambda1/param.m r, 0 0, 0 0, -param.lambda1/param.m*[cos(psit) sin(psit)] ;...
    -param.g(3) 0, -r -param.lambda1/param.m, 0 0, 0 0, -param.lambda1/param.m*[-sin(psit) cos(psit)];...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];
% x = [phi theta u v x y biasphi biastheta]
F6 = @(psit, r) eye(8)+param.dt*[  0 r, 0 0, 0 0, -1 0;...
    -r 0, 0 0, 0 0, 0 -1;...
    0 param.g(3), 0 r, 0 0, 0 0;...
    -param.g(3) 0, -r 0, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];

B5 = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0]';
B6 = param.dt*[1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0 -1 0 0 0 0 0;...
    0 0 0 -1 0 0 0 0]';

G5 = [1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0;...
    0 0, param.dt^2/2 0, 0 0, 0 0;...
    0 0, 0 param.dt^2/2, 0 0, 0 0;...
    0 0, 0 0, 1 0 0 0;...
    0 0, 0 0, 0 1 0 0;...
    0 0, 0 0, 0 0 1 0;...
    0 0, 0 0, 0 0 0 1];
G6 = [1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0;...
    0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 1];

H51 = @(psit) [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0,cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0];

H52 = @(psit) [0 0, param.lambda1/param.m 0, 0 0, 0 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
    0 0, 0 param.lambda1/param.m, 0 0, 0 0, param.lambda1/param.m*[-sin(psit) cos(psit)]];

R51 = diag([Rpos Rpos Rspeed Rspeed]);...
    R52 = diag([Racc Racc]);
R5 = blkdiag(R51,R52);
R6 = R51;
Q5 = diag([Qeul, Qv, Qbiasg, Qwind]);
Q6 = diag([Qeul, Qv, Qx, Qbiasg]);
% Stationary Kalman

x5 = [zeros(8,1); 0*windData(1,1:2)']; x6 = zeros(8,1);
P5 = 1000*eye(10);P6 = 1000*eye(8);
for ii = 1:length(clock)
    % Data
    u5 = [pnoise(ii)+bias_p(ii) qnoise(ii)+bias_q(ii)];
    y5 = [xnoise(ii) ynoise(ii),...
        xdotnoise(ii), ydotnoise(ii),...
        axnoise(ii) aynoise(ii) ];
    u6 = [u5 axnoise(ii) aynoise(ii)];
    y6 = [xnoise(ii) ynoise(ii),...
        xdotnoise(ii), ydotnoise(ii)];
    
    % Prediction step
    F5ii = F5(psinoise(ii),psidot(ii));
    x5 = F5ii*x5 + B5*u5';
    P5 = F5ii*P5*F5ii' + G5*Q5*G5';
    
    F6ii = F6(psinoise(ii),psidot(ii));
    x6 = F6ii*x6 + B6*u6';
    P6 = F6ii*P6*F6ii' + G6*Q6*G6';
    
    H51ii =  H51(psinoise(ii));
    H52ii = H52(psinoise(ii));
    H5 = [H51ii;H52ii];
    H6 = H51ii(:,1:8);
    
    rr5(ii) = rank(obsv(F5ii,H5));
    rr6(ii) = rank(obsv(F6ii,H6));
    cond5(ii) = cond(obsv(F5ii,H5));
    cond6(ii) = cond(obsv(F6ii,H6));
    
    if mod(ii-1,index) == 0
        % Correction step
        S5 = H5*P5*H5'+R5;
        K5 = P5*H5'*inv(S5);
        x5 = x5 + K5*(y5' - H5*x5);
        P5 = P5 - K5*H5*P5;
        
        % Correction step
        S6 = H6*P6*H6'+R6;
        K6 = P6*H6'*inv(S6);
        x6 = x6 + K6*(y6' - H6*x6);
        P6 = P6 - K6*H6*P6;
    else
        % Correction step
        S5 =H52ii*P5*H52ii'+R52;
        K5 = P5*H52ii'*inv(S5);
        x5 = x5 + K5*(y5(5:end)' - H52ii*x5);
        P5 = P5 - K5*H52ii*P5;
        
    end
    
    P5m(ii,:) = diag(P5);
    P6m(ii,:) = diag(P6);
    x5m(ii,:) = x5;
    y5m(ii,:) = H5*x5;
    x6m(ii,:) = x6;
    y6m(ii,:) = H6*x6;
end
disp('Done compute Gyro bias')

%%
eulLim = 0.5;beulLim = 0.5;windxLim = 2;windyLim = 2;
horix = -1900 ;very = 700;movex = 450; movey = 450; figx = 420; figy = 320;%560 420
close all
ff(1) = figure(1);clf
ff(1).Position = [horix very figx figy];
subplot(2,1,1)
plot(clock,phinoise,'linewidth',1.5)
hold on
plot(clock,x5m(:,1),'-','linewidth',1.5)
plot(clock,x6m(:,1),'-.','linewidth',1.5,'color','g')
plot(clock,0.1*psinoise,'linewidth',1.5)
ylim([-eulLim eulLim])

legend 'phi true' 'phi Gyro+aout' 'phi Gyro+ain'% 'theta est p+a'
subplot(2,1,2)
plot(clock,thetanoise,'linewidth',1.5)
hold on
plot(clock,x5m(:,2),'-','linewidth',1.5)
plot(clock,x6m(:,2),'-.','linewidth',1.5,'color','g')
plot(clock,0.1*psinoise,'linewidth',1.5)
legend 'theta true' 'theta Gyro+aout' 'theta Gyro+ain'%'theta est p+a'
ylim([-eulLim eulLim])

ff(2) = figure(2);clf
ff(2).Position = [horix+movex very figx figy];
subplot(2,1,1)
plot(clock,bias_p,'linewidth',1.5)
hold on
plot(clock,x5m(:,7),'-','linewidth',1.5)
plot(clock,x6m(:,7),'-.','linewidth',1.5,'color','g')
legend 'bias_p' 'Gyro+aout' 'Gyro ain'
ylim([-beulLim beulLim])

subplot(2,1,2)
plot(clock,bias_q,'linewidth',1.5)
hold on
plot(clock,x5m(:,8),'-','linewidth',1.5)
plot(clock,x6m(:,8),'-.','linewidth',1.5,'color','g')
legend 'bias_q' 'Gyro+aout' 'Gyro ain'
ylim([-beulLim beulLim])


%
ff(3) = figure(3);clf
ff(3).Position = [horix+2*movex very figx figy];
subplot(2,1,1)
plot(clock,windData(:,1),'linewidth',1.5)
hold on
plot(clock,x5m(:,9),'-','linewidth',1.5)
legend 'u_w' 'u_w Gyro+aout'
ylim([-windxLim windxLim])

subplot(2,1,2)
plot(clock,windData(:,2),'linewidth',1.5)
hold on
plot(clock,x5m(:,10),'-','linewidth',1.5)
legend 'v_w' 'v_w Gyro+aout'
ylim([-windyLim windyLim])


ff(4) = figure(4);clf
ff(4).Position = [horix+3*movex very figx figy];
subplot(2,1,1)
plot(clock,unoise,'linewidth',1.5)
hold on
plot(clock,x5m(:,3),'-','linewidth',1.5)
plot(clock,x6m(:,3),'-.','linewidth',1.5,'color','g')
legend 'xdot' 'xdot Gyro+aout' 'xdot Gyro+ain'
ylim([-windxLim windxLim])

subplot(2,1,2)
plot(clock,vnoise,'linewidth',1.5)
hold on
plot(clock,x5m(:,4),'-','linewidth',1.5)
plot(clock,x6m(:,4),'-.','linewidth',1.5,'color','g')
legend 'ydot' 'ydot Gyro+aout' 'ydot Gyro+ain'
ylim([-windyLim windyLim])

%%
% Estimate acceleration bias without euler angle measurement, we use yaw
% angle change instead of
% x = [phi theta u v x y biasx biasy windx windy]
index = 200
% Tuning coefficient
Qeul = 1e-10*[1 1];
Qv = 1e-10*[1 1];
Qx = 1e-10*[1 1];
Qbiasg = 5*1e-5*[1 1];
Qbiasa = 5*1e-6*[1 1];
Qwind = 1e-4*[1 1];
index = 200

close all
clear x7m y7m x8m y8m
F7 = @(psit, r) eye(10)+param.dt*[  0 r, 0 0, 0 0, 0 0, 0 0;...
    -r 0, 0 0, 0 0, 0 0, 0 0;...
    0 param.g(3), -param.lambda1/param.m r, 0 0, 0 0, -param.lambda1/param.m*[cos(psit) sin(psit)] ;...
    -param.g(3) 0, -r -param.lambda1/param.m, 0 0, 0 0, -param.lambda1/param.m*[-sin(psit) cos(psit)];...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];

% x = [phi theta u v x y biasx biasy]
F8 = @(psit, r) eye(8)+param.dt*[  0 r, 0 0, 0 0, 0 0;...
    -r 0, 0 0, 0 0, 0 0;...
    0 param.g(3), 0 r, 0 0, 1 0;...
    -param.g(3) 0, -r 0, 0 0, 0 1;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H71 = @(psit) [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0];
H72 = @(psit) [0 0, param.lambda1/param.m 0, 0 0, 1 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
    0 0, 0 param.lambda1/param.m, 0 0, 0 1, param.lambda1/param.m*[-sin(psit) cos(psit)]];

H8 = @(psit) [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0];

B7 = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0]';
B8 = param.dt*[1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, -1 0, 0 0, 0 0;...
    0 0, 0 -1, 0 0, 0 0]';
G7 = [1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0, 0 0;...
    0 0, 0 0, 0 0, 1 0 0 0;...
    0 0, 0 0, 0 0, 0 1 0 0;...
    0 0, 0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 0, 0 1];
G8 = [1 0, 0 0,0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0;...
    0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 1];

R71 = diag([Rpos Rpos  Rspeed Rspeed]);...
    R72 = diag([Racc Racc]);
R7 = blkdiag(R71,R72);
R8 = R71;
Q7 = diag([Qeul, Qv, Qx, Qbiasa, Qwind]);
Q8 = diag([Qeul, Qv, Qx, Qbiasa]);
% Stationary Kalman
% [K7bar,Pp7,Pf7] = dlqe(F7,G7,H7,Q7,R7);

x7 = [zeros(8,1);0*windData(1,1:2)']; x8 = zeros(8,1);
P7 = 1000*eye(10); P8 = 1000*eye(8);
for ii = 1:length(clock)
    % Data
    u7 = [pnoise(ii)+0*bias_p(ii) qnoise(ii)+0*bias_q(ii)];
    y7 = [xnoise(ii) ynoise(ii),...
        xdotnoise(ii) ydotnoise(ii),...
        axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii) ];
    u8 = [pnoise(ii)+0*bias_p(ii) qnoise(ii)+0*bias_q(ii) axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii) ];
    y8 = [xnoise(ii) ynoise(ii),...
        xdotnoise(ii) ydotnoise(ii)];
    
    F7ii = F7(psinoise(ii),psidot(ii));
    H71ii = H71(psinoise(ii));
    H72ii = H72(psinoise(ii));
    H7 = [H71ii;H72ii];
    
    F8ii = F8(psinoise(ii),psidot(ii));
    H8ii = H8(psinoise(ii));
    
    
    % Prediction step
    x7 = F7ii*x7 + B7*u7';
    P7 = F7ii*P7*F7ii' + G7*Q7*G7';
    
    x8 = F8ii*x8 + B8*u8';
    P8 = F8ii*P8*F8ii' + G8*Q8*G8';
    
    rr7(ii) = rank(obsv(F7ii,H7));
    rr8(ii) = rank(obsv(F8ii,H8ii));
    cond7(ii) = cond(obsv(F7ii,H7));
    cond8(ii) = cond(obsv(F8ii,H8ii));
    
    if mod(ii-1,index) == 0
        % Correction step
        S7 = H7*P7*H7'+R7;
        K7 = P7*H7'*inv(S7);
        x7 = x7 + K7*(y7' - H7*x7);
        P7 = P7 - K7*H7*P7;
        
        S8 = H8ii*P8*H8ii'+R8;
        K8 = P8*H8ii'*inv(S8);
        x8 = x8 + K8*(y8' - H8ii*x8);
        P8 = P8 - K8*H8ii*P8;
    else
        % Correction step
        S7 = H72ii*P7*H72ii'+R72;
        K7 = P7*H72ii'*inv(S7);
        x7 = x7 + K7*(y7(5:end)' - H72ii*x7);
        P7 = P7 - K7*H72ii*P7;
        
    end
    
    
    x7m(ii,:) = x7;
    y7m(ii,:) = H7*x7;
    
    x8m(ii,:) = x8;
    y8m(ii,:) = H8ii*x8;
end
disp('Done compute acc bias ')


%
baccLim = 0.6;windxLim = 2;windyLim = 2;uLim = 3;vLim = 3;

ff(1) = figure(1);clf
ff(1).Position = [horix very figx figy];
subplot(2,1,1)
plot(clock,phinoise,'linewidth',1.5)
hold on
plot(clock,x7m(:,1),'-','linewidth',1.5)
plot(clock,x8m(:,1),'-.','linewidth',1.5,'color','g')
plot(clock,0.1*psinoise,'linewidth',1.5)
ylim([-eulLim eulLim])

legend 'phi true' 'phi Gyro' 'phi Gyro+Non'% 'theta est p+a'
subplot(2,1,2)
plot(clock,thetanoise,'linewidth',1.5)
hold on
plot(clock,x7m(:,2),'-','linewidth',1.5)
plot(clock,x8m(:,2),'-.','linewidth',1.5,'color','g')
plot(clock,0.1*psinoise,'linewidth',1.5)
legend 'theta true' 'theta Gyro' 'theta Gyro+Non'%'theta est p+a'
ylim([-eulLim eulLim])

ff(2) = figure(2);clf
ff(2).Position = [horix+movex very figx figy];
subplot(2,1,1)
plot(clock,unoise,'linewidth',1.5)
hold on
plot(clock,x7m(:,3),'-','linewidth',1.5)
plot(clock,x8m(:,3),'-.','linewidth',1.5,'color','g')
legend 'bias ax true' 'bias ax aout' 'bias ax aout+non'
ylim([-uLim uLim])

subplot(2,1,2)
plot(clock,vnoise,'linewidth',1.5)
hold on
plot(clock,x7m(:,4),'-','linewidth',1.5)
plot(clock,x8m(:,4),'-.','linewidth',1.5,'color','g')
legend 'bias ay true'  'bias ay aout' 'bias ay aout+non'
ylim([-vLim vLim])


ff(4) = figure(4);clf
ff(4).Position = [horix+3*movex very figx figy];
subplot(2,1,1)
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,x7m(:,7),'-','linewidth',1.5)
plot(clock,x8m(:,7),'-.','linewidth',1.5,'color','g')
legend 'bias ax true' 'bias ax aout' 'bias ax aout+non'
ylim([-baccLim baccLim])

subplot(2,1,2)
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,x7m(:,8),'-','linewidth',1.5)
plot(clock,x8m(:,8),'-.','linewidth',1.5,'color','g')
legend 'bias ay true'  'bias ay aout' 'bias ay aout+non'
ylim([-baccLim baccLim])

ff(5) = figure(5);clf
ff(5).Position = [horix very-movey figx figy];
subplot(2,1,1)
plot(clock,windData(:,1),'linewidth',1.5)
hold on
plot(clock,x7m(:,9),'-','linewidth',1.5)
legend 'windx true' 'windx aout'
ylim([-windxLim windxLim])

subplot(2,1,2)
plot(clock,windData(:,2),'linewidth',1.5)
hold on
plot(clock,x7m(:,10),'-','linewidth',1.5)
legend 'windy true'  'windy aout'
ylim([-windyLim windyLim])


ff(6) = figure(6);clf
ff(6).Position = [horix+movex very-movey figx figy];
subplot(2,1,1)
plot(clock,abs(stateData(:,4)-x7m(:,1)),'linewidth',1.5)
hold on
plot(clock,abs(stateData(:,4)-x8m(:,1)),'linewidth',1.5)
legend 'Diff phi acc aout' 'Diff phi acc aout+non'
ylim([-eulLim eulLim])
subplot(2,1,2)
plot(clock,abs(stateData(:,5)-x7m(:,2)),'linewidth',1.5)
hold on
plot(clock,abs(stateData(:,5)-x8m(:,2)),'linewidth',1.5)
legend 'Diff theta acc aout' 'Diff theta acc aout+non'
ylim([-eulLim eulLim])


%% Now we apply multiple model fault detection in compensation for mass change detection
% WE need to defined three model which is fault free, gyro fault and acc
% fault

index = 200;
% Tuning coefficient
% Try to tune for the bias and then tune wind variance
Qeul = 1e-3*1e0*1e-5*[1 1]; % 1e-7 make no detection, 1e-9 false alarm when acc bias,
Qv = 1e-3*1e7*1e-5*[1 1];
Qx = 1e-3*1e7*1e-5*[1 1];
Qbiasg = 1e-3*1e2*1e-4*[1 1]; % 10-10 10-10 10-10 5 10-5 5 10-5 10-4
Qbiasa = 1e-3*1e2*5*1e-5*[1 1];
Qwind = 5*1e-3*[1 1];

% Good combo
% 1e-8 1e-8 1e-8  51e-5 1e-5 1e-4
% 1e-10 1e-10 1e-10 51e-5 51e-5 1e-4
% Design yaw rate is important, yaw should change fast enough to excite some
% dynamic in the system
% Adding the gain to change noise bias levels
gain_bp = 0.0*ones(N,1); % 1-10 -10 -10 5-6 5-5 -4
gain_bq = 0.0*ones(N,1);
gain_bx = 1.0*ones(N,1); % 5 10-5 5 10-5 10-4
gain_by = 1.0*ones(N,1); % 5 1e-5 1e-5 1e-4

% Qbiasg = 5*1e-7*[1 1];
% Qbiasa = 5*1e-6*[1 1];
% Qwind = 1e-4*[1 1];



% Fault free model x = [phi theta u v x y windx windy];
F0out = @(psit, r) eye(8)+param.dt*[0 r, 0 0, 0 0, 0 0;...
    -r 0, 0 0, 0 0, 0 0;...
    0 param.g(3), -param.lambda1/param.m r, 0 0, -param.lambda1/param.m*[cos(psit) sin(psit)] ;...
    -param.g(3) 0, -r -param.lambda1/param.m, 0 0, -param.lambda1/param.m*[-sin(psit) cos(psit)];...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
B0out = param.dt*[1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0]';

% x = [phi theta u v x y]
F0in = @(psit, r) eye(6)+param.dt*[  0 r, 0 0, 0 0;...
    -r 0, 0 0, 0 0;...
    0 param.g(3), 0 r, 0 0;...
    -param.g(3) 0, -r 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0;...
    0 0, sin(psit) cos(psit), 0 0];
B0in = param.dt*[1 0, 0 0, 0 0;...
    0 1, 0 0, 0 0;...
    0 0, -1 0, 0 0;...
    0 0, 0 -1, 0 0]';
G0out = [1 0, 0 0, 0 0, 0 0 ;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0;...
    0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 1];
G0in = [1 0, 0 0, 0 0;...
    0 1, 0 0, 0 0;...
    0 0, param.dt 0, 0 0;...
    0 0, 0 param.dt, 0 0;...
    0 0, 0 0, param.dt^2/2 0;...
    0 0, 0 0, 0 param.dt^2/2];
% GPS
H01out = @(psit) [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0];
% ACC
H02out = @(psit) [0 0, param.lambda1/param.m 0, 0 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
    0 0, 0 param.lambda1/param.m, 0 0, param.lambda1/param.m*[-sin(psit) cos(psit)]];

H0in = @(psit) [0 0, 0 0, 1 0;...
    0 0, 0 0, 0 1;...
    0 0, cos(psit) -sin(psit), 0 0;...
    0 0, sin(psit) cos(psit), 0 0];

R01out = diag([Rpos Rpos Rspeed Rspeed]);...
    R02out = diag([Racc Racc]);
R0out = blkdiag(R01out,R02out);
R0in = R01out;
Q0out = diag([Qeul, Qv, Qx, Qwind]);
Q0in = diag([Qeul, Qv, Qx]);



% Gyro fault model
% x = [phi theta u v x y biasphi biastheta windx windy]
Fgout = @(psit, r) eye(10)+param.dt*[0 r, 0 0, 0 0, -1 0, 0 0;...
    -r 0, 0 0, 0 0, 0 -1, 0 0;...
    0 param.g(3), -param.lambda1/param.m r, 0 0, 0 0, -param.lambda1/param.m*[cos(psit) sin(psit)] ;...
    -param.g(3) 0, -r -param.lambda1/param.m, 0 0, 0 0, -param.lambda1/param.m*[-sin(psit) cos(psit)];...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];
% x = [phi theta u v x y biasphi biastheta]
Fgin = @(psit, r) eye(8)+param.dt*[  0 r, 0 0, 0 0, -1 0;...
    -r 0, 0 0, 0 0, 0 -1;...
    0 param.g(3), 0 r, 0 0, 0 0;...
    -param.g(3) 0, -r 0, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];

Bgout = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0]';
Bgin = param.dt*[1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0 -1 0 0 0 0 0;...
    0 0 0 -1 0 0 0 0]';

Ggout = [1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0, 0 0;...
    0 0, 0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 0, 0 1, 0 0;...
    0 0, 0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 0, 0 1];
Ggin = [1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0;...
    0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 1];

% GPS
Hg1out = @(psit) [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0];
% ACC
Hg2out = @(psit) [0 0, param.lambda1/param.m 0, 0 0, 0 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
    0 0, 0 param.lambda1/param.m, 0 0, 0 0, param.lambda1/param.m*[-sin(psit) cos(psit)]];

Hgin = @(psit) [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0];

Rg1out = diag([Rpos Rpos Rspeed Rspeed]);...
    Rg2out = diag([Racc Racc]);
Rgout = blkdiag(Rg1out,Rg2out);
Rgin = Rg1out;
Qgout = diag([Qeul, Qv, Qx, Qbiasg, Qwind]);
Qgin = diag([Qeul, Qv, Qx, Qbiasg]);

% ACC fault model

Faout = @(psit, r) eye(10)+param.dt*[  0 r, 0 0, 0 0, 0 0, 0 0;...
    -r 0, 0 0, 0 0, 0 0, 0 0;...
    0 param.g(3), -param.lambda1/param.m r, 0 0, 0 0, -param.lambda1/param.m*[cos(psit) sin(psit)] ;...
    -param.g(3) 0, -r -param.lambda1/param.m, 0 0, 0 0, -param.lambda1/param.m*[-sin(psit) cos(psit)];...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];

% x = [phi theta u v x y biasx biasy]
Fain = @(psit, r) eye(8)+param.dt*[  0 r, 0 0, 0 0, 0 0;...
    -r 0, 0 0, 0 0, 0 0;...
    0 param.g(3), 0 r, 0 0, 1 0;...
    -param.g(3) 0, -r 0, 0 0, 0 1;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
Ha1out = @(psit) [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0];
Ha2out = @(psit) [0 0, param.lambda1/param.m 0, 0 0, 1 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
    0 0, 0 param.lambda1/param.m, 0 0, 0 1, param.lambda1/param.m*[-sin(psit) cos(psit)]];

Hain = @(psit) [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0];

Baout = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0]';
Bain = param.dt*[1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, -1 0, 0 0, 0 0;...
    0 0, 0 -1, 0 0, 0 0]';
Gaout = [1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0, 0 0;...
    0 0, 0 0, 0 0, 1 0 0 0;...
    0 0, 0 0, 0 0, 0 1 0 0;...
    0 0, 0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 0, 0 1];
Gain = [1 0, 0 0,0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0;...
    0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 1];

Ra1out = diag([Rpos Rpos  Rspeed Rspeed]);...
    Ra2out = diag([Racc Racc]);
Raout = blkdiag(Ra1out,Ra2out);
Rain = Ra1out;
Qaout = diag([Qeul, Qv, Qx, Qbiasa, Qwind]);
Qain = diag([Qeul, Qv, Qx, Qbiasa]);
% Stationary Kalman
% [K7bar,Pp7,Pf7] = dlqe(F7,G7,H7,Q7,R7);

%
Rpos = 9*0.01; Rspeed = 0.01; Reul = 0.01; Racc  = 1e-6; Rr = 1e-4;Rpsi = 5*1e-5;
close all
clear x7m y7m x8m y8m



Nmc = 4; % We consider r^{mc}
tic
for mc = 1:Nmc
    
    % Stress the mothod based on psinoise
    xnoise = stateData(:,1) + sqrt(Rpos)*randn(N,1);
    ynoise = stateData(:,2) + sqrt(Rpos)*randn(N,1);
    xdotnoise = stateData(:,7) + sqrt(Rspeed)*randn(N,1);
    ydotnoise = stateData(:,8) + sqrt(Rspeed)*randn(N,1);
    
    Rr = 10^(-(3*mc-2));Rpsi = 5*10^(-(3*mc-1));
    % Rr = 1e-4;Rpsi = 5*1e-5;
    psinoise = stateData(1:end,6) + sqrt(Rpsi)*randn(N,1);
    psidot = stateData(1:end,18) + sqrt(Rr)*randn(N,1); % Similar performance to state 18, r
    xaout = [zeros(8,1);0;0]; xain = zeros(8,1);
    Paout = 1000*eye(10); Pain = 1000*eye(8);
    
    
    % Start Kalman filter
    % a as output
    x0out = zeros(8,1);P0out = 1000*eye(8);
    xgout = zeros(10,1);Pgout = 1000*eye(10);
    % xaout = zeros(10,1);Paout = 1000*eye(10);
    
    % a as input
    x0in = zeros(6,1);P0in = 1000*eye(6);
    xgin = zeros(8,1);Pgin = 1000*eye(8);
    % xain = zeros(8,1);Pain = 1000*eye(8);
    
    P_0out = 0.99;P_gout = 0.005; P_aout = 0.005;
    P_0in = 0.99;P_gin = 0.005; P_ain = 0.005;
    P_0out_old = 0.99;P_gout_old = 0.005; P_aout_old = 0.005;
    P_0in_old = 0.99;P_gin_old = 0.005; P_ain_old = 0.005;
    w0out = 0; wgout = 0; waout = 0;
    w0in = 0; wgin = 0; wain = 0;
    P_thres = 1e-8;
    
    
    for ii = 1:length(clock)
        % Preperation step
        uout = [pnoise(ii)+gain_bp(ii)*bias_p(ii) qnoise(ii)+gain_bq(ii)*bias_q(ii)]';
        yout = [xnoise(ii) ynoise(ii) xdotnoise(ii) ydotnoise(ii),...
            axnoise(ii)+gain_bx(ii)*bias_ax(ii) aynoise(ii)+gain_by(ii)*bias_ay(ii)]';
        
        uin = [pnoise(ii)+gain_bp(ii)*bias_p(ii) qnoise(ii)+gain_bq(ii)*bias_q(ii),...
            axnoise(ii)+gain_bx(ii)*bias_ax(ii) aynoise(ii)+gain_by(ii)*bias_ay(ii)]';
        yin = [xnoise(ii) ynoise(ii),...
            xdotnoise(ii) ydotnoise(ii)]';
        
        
        % Prediction step
        F0outii = F0out(psinoise(ii),psidot(ii));
        Fgoutii = Fgout(psinoise(ii),psidot(ii));
        Faoutii = Faout(psinoise(ii),psidot(ii));
        
        F0inii = F0in(psinoise(ii),psidot(ii));
        Fginii = Fgin(psinoise(ii),psidot(ii));
        Fainii = Fain(psinoise(ii),psidot(ii));
        
        % Update state
        x0out = F0outii*x0out + B0out*uout;
        xgout = Fgoutii*xgout + Bgout*uout;
        xaout = Faoutii*xaout + Baout*uout;
        
        x0in = F0inii*x0in + B0in*uin;
        xgin = Fginii*xgin + Bgin*uin;
        xain = Fainii*xain + Bain*uin;
        
        % Update covariance
        P0out = F0outii*P0out*F0outii' + G0out*Q0out*G0out';
        Pgout = Fgoutii*Pgout*Fgoutii' + Ggout*Qgout*Ggout';
        Paout = Faoutii*Paout*Faoutii' + Gaout*Qaout*Gaout';
        
        P0in = F0inii*P0in*F0inii' + G0in*Q0in*G0in';
        Pgin = Fginii*Pgin*Fginii' + Ggin*Qgin*Ggin';
        Pain = Fainii*Pain*Fainii' + Gain*Qain*Gain';
        
        
        
        % Measurement update step
        H01outii = H01out(psinoise(ii));
        H02outii = H02out(psinoise(ii));
        H0outii = [H01outii; H02outii];
        H0inii = H0in(psinoise(ii));
        
        Hg1outii = Hg1out(psinoise(ii));
        Hg2outii = Hg2out(psinoise(ii));
        Hgoutii = [Hg1outii; Hg2outii];
        Hginii = Hgin(psinoise(ii));
        
        Ha1outii = Ha1out(psinoise(ii));
        Ha2outii = Ha2out(psinoise(ii));
        Haoutii = [Ha1outii;Ha2outii];
        Hainii = Hain(psinoise(ii));
        
        
        
        if mod(ii-1,index) == 0 % Update when GPS availabel
            
            % Fault free
            S0out = H0outii*P0out*H0outii'+R0out;
            invS0out = eye(size(S0out))/S0out;
            K0out = P0out*H0outii'*invS0out;
            r0out = yout - H0outii*x0out;
            x0out = x0out + K0out*r0out;
            P0out = P0out - K0out*H0outii*P0out;
            
            S0in = H0inii*P0in*H0inii'+R0in;
            invS0in = eye(size(S0in))/S0in;
            K0in = P0in*H0inii'*invS0in;
            r0in = yin - H0inii*x0in;
            x0in = x0in + K0in*r0in;
            P0in = P0in - K0in*H0inii*P0in;
            
            % Gyro fault
            Sgout = Hgoutii*Pgout*Hgoutii'+Rgout;
            invSgout = eye(size(Sgout))/Sgout;
            Kgout = Pgout*Hgoutii'*invSgout;
            rgout = yout - Hgoutii*xgout;
            xgout = xgout + Kgout*rgout;
            Pgout = Pgout - Kgout*Hgoutii*Pgout;
            
            Sgin = Hginii*Pgin*Hginii'+Rgin;
            invSgin = eye(size(Sgin))/Sgin;
            Kgin = Pgin*Hginii'*invSgin;
            rgin = yin - Hginii*xgin;
            xgin = xgin + Kgin*rgin;
            Pgin = Pgin - Kgin*Hginii*Pgin;
            
            % ACC fault
            Saout = Haoutii*Paout*Haoutii'+Raout;
            invSaout = eye(size(Saout))/Saout;
            Kaout = Paout*Haoutii'*invSaout;
            raout = yout - Haoutii*xaout;
            xaout = xaout + Kaout*raout;
            Paout = Paout - Kaout*Haoutii*Paout;
            
            Sain = Hainii*Pain*Hainii'+Rain;
            invSain = eye(size(Sain))/Sain;
            Kain = Pain*Hainii'*invSain;
            rain = yin - Hainii*xain;
            xain = xain + Kain*rain;
            Pain = Pain - Kain*Hainii*Pain;
            
            
            % Compute the weight of the wieght
            w0out = 1/sqrt(det(S0out))*exp(-0.5*r0out'*invS0out*r0out);
            wgout = 1/sqrt(det(Sgout))*exp(-0.5*rgout'*invSgout*rgout);
            waout = 1/sqrt(det(Saout))*exp(-0.5*raout'*invSaout*raout);
            w0in = 1/sqrt(det(S0in))*exp(-0.5*r0in'*invS0in*r0in);
            wgin = 1/sqrt(det(Sgin))*exp(-0.5*rgin'*invSgin*rgin);
            wain = 1/sqrt(det(Sain))*exp(-0.5*rain'*invSain*rain);
            
            % Store the previous weight
            P_0out_old = P_0out;
            P_gout_old = P_gout;
            P_aout_old = P_aout;
            P_0in_old = P_0in;
            P_gin_old = P_gin;
            P_ain_old = P_ain;
            
            % Update the weight
            all_weight_out = w0out*P_0out_old + wgout*P_gout_old+waout*P_aout_old;
            P_0out = w0out*P_0out_old/all_weight_out;
            P_gout = wgout*P_gout_old/all_weight_out;
            P_aout = waout*P_aout_old/all_weight_out;
            
            all_weight_in = w0in*P_0in_old + wgin*P_gin_old + wain*P_ain_old;
            P_0in = w0in*P_0in_old/all_weight_in;
            P_gin = wgin*P_gin_old/all_weight_in;
            P_ain = wain*P_ain_old/all_weight_in;
            
            %         if P_0out<P_thres
            %             P_0out = P_thres;
            %         end
            %         if P_gout<P_thres
            %             P_gout = P_thres;
            %         end
            %         if P_aout<P_thres
            %             P_aout = P_thres;
            %         end
            %
            %         if P_0in<P_thres
            %             P_0in = P_thres;
            %         end
            %         if P_gin<P_thres
            %             P_gin = P_thres;
            %         end
            %         if P_ain<P_thres
            %             P_ain = P_thres;
            %         end
            
            
            
            %         % Compute the rank of the matrix
            %         % Estimate the observability matrix
            %         if ii<N-200
            %             [~,rank0out] = quadObsv(F0out,H01out,H02out,psinoise(ii-199:ii),psidot(ii-199:ii));
            %             [~,rankgout] = quadObsv(Fgout,Hg1out,Hg2out,psinoise(ii-199:ii),psidot(ii-199:ii));
            %             [~,rankaout] = quadObsv(Faout,Ha1out,Ha2out,psinoise(ii-199:ii),psidot(ii-199:ii));
            %             [~,rank0in] = quadObsv(F0in,H0in,psinoise(ii-199:ii),psidot(ii-199:ii));
            %             [~,rankgin] = quadObsv(Fgin,Hgin,psinoise(ii-199:ii),psidot(ii-199:ii));
            %             [~,rankain] = quadObsv(Fain,Hain,psinoise(ii-199:ii),psidot(ii-199:ii));
            %         end
            %         rr0out(ii,:) = rank0out;
            %         rrgout(ii,:) = rankgout;
            %         rraout(ii,:) = rankaout;
            %         rr0in(ii,:) = rank0in;
            %         rrgin(ii,:) = rankgin;
            %         rrain(ii,:) = rankain;
            
        else % Update at 200Hz without GPS measurement, ACC only
            % Fault free
            S0out2 = H02outii*P0out*H02outii'+R02out;
            invS0out2 = S0out2\eye(size(S0out2));
            K0out = P0out*H02outii'*invS0out2;
            r0out2 = yout(5:end) - H02outii*x0out;
            x0out = x0out + K0out*r0out2;
            P0out = P0out - K0out*H02outii*P0out;
            
            %         K0in = P0in*H02inii'*inv(H02inii*P0in*H02inii'+R0in);
            %         x0in = x0in + K0in*(yin(5:end) - H02inii*x0in);
            %         P0in = P0in - K0in*H02inii*P0in;
            
            % Gyro fault
            Sgout2 = Hg2outii*Pgout*Hg2outii'+Rg2out;
            invSgout2 = Sgout2\eye(size(Sgout2));
            Kgout = Pgout*Hg2outii'*invSgout2;
            rgout2 = yout(5:end) - Hg2outii*xgout;
            xgout = xgout + Kgout*rgout2;
            Pgout = Pgout - Kgout*Hg2outii*Pgout;
            
            %         Kgin = Pgin*Hg2inii'*inv(Hginii*Pgin*Hginii'+Rgin);
            %         xgin = xgin + Kgin*(yin(5:end) - Hginii*xgin);
            %         Pgin = Pgin - Kgin*Hginii*Pgin;
            
            % ACC fault
            Saout2 = Ha2outii*Paout*Ha2outii'+Ra2out;
            invSaout2 = Saout2\eye(size(Saout2));
            Kaout = Paout*Ha2outii'*invSaout2;
            raout2 = yout(5:end) - Ha2outii*xaout;
            xaout = xaout + Kaout*raout2;
            Paout = Paout - Kaout*Ha2outii*Paout;
            
            %         Kain = Pain*Hainii'*inv(Hainii*Pain*Hainii'+Rain);
            %         xain = xain + Kain*(yin(5:end) - Hainii*xain);
            %         Pain = Pain - Kain*Hainii*Pain;
            %
            
        end
        
        
        
        
        x0outm(ii,:,mc) = x0out;
        y0outm(ii,:,mc) = H0outii*x0out;
        
        x0inm(ii,:,mc) = x0in;
        y0inm(ii,:,mc) = H0inii*x0in;
        
        xgoutm(ii,:,mc) = xgout;
        ygoutm(ii,:,mc) = Hgoutii*xgout;
        
        xginm(ii,:,mc) = xgin;
        yginm(ii,:,mc) = Hginii*xgin;
        
        xaoutm(ii,:,mc) = xaout;
        yaoutm(ii,:,mc) = Haoutii*xaout;
        
        xainm(ii,:,mc) = xain;
        yainm(ii,:,mc) = Hainii*xain;
        
        % Store probability
        P_0outm(ii,:,mc) = P_0out;
        P_goutm(ii,:,mc) = P_gout;
        P_aoutm(ii,:,mc) = P_aout;
        
        w0outm(ii,:,mc) = w0out;
        wgoutm(ii,:,mc) = wgout;
        waoutm(ii,:,mc) = waout;
        w0inm(ii,:,mc) = w0in;
        wginm(ii,:,mc) = wgin;
        wainm(ii,:,mc) = wain;
        
        P_0inm(ii,:,mc) = P_0in;
        P_ginm(ii,:,mc) = P_gin;
        P_ainm(ii,:,mc) = P_ain;
        
        xoutm(ii,:,mc) = P_0out*x0out(1:6) + P_gout*xgout(1:6) + P_aout*xaout(1:6);
        xinm(ii,:,mc) = P_0in*x0in(1:6) + P_gin*xgin(1:6) + P_ain*xain(1:6);
        bias_gyroout(ii,:,mc) =  P_gout*xgout(7:8);
        bias_accout(ii,:,mc) =  P_aout*xaout(7:8);
        
        bias_gyroin(ii,:,mc) =  P_gin*xgin(7:8);
        bias_accin(ii,:,mc) =  P_ain*xain(7:8);
        
    end
    
    mc
    toc
    tic
end



disp('Done compute multiple adaptive fault detection')

toc
%
eulLim = 0.75;beulLim = 0.75;vxLim = 2;vyLim = 2;windxLim = 4;windyLim = 2;
horix = -1900 ;very = 700;movex = 450; movey = 450; figx = 420; figy = 320;%560 420
% -1950 700 420 320 450 300
close all
if Nmc == 1
    xmc = 1;
else
    xmc = randi(Nmc,1,1);
end

ff(1) = figure(1);clf
ff(1).Position = [horix very-0.5*movey 1.5*figx 1.5*figy];
h2 = axes;h2.FontSize = 14;hold on
subplot(2,1,1);
plot(clock,(phinoise-squeeze(xoutm(:,1,xmc))),'linewidth',1.5)
hold on
plot(clock,(phinoise-squeeze(xinm(:,1,xmc))),'-.','linewidth',1.5)
legend({'$\phi - \hat{\phi}$ Ayaw MMAE', '$\phi - \hat{\phi}$ Byaw MMAE'},'interpreter','latex','FontSize',14)
ylabel('$|\phi-\hat{\phi}|$','interpreter','latex','FontSize',14)
xlabel('Time [s]','FontSize',14)
ylim([-beulLim beulLim])
grid on

subplot(2,1,2)
plot(clock,(thetanoise-squeeze(xoutm(:,2,xmc))),'-','linewidth',1.5)
hold on
plot(clock,(thetanoise-squeeze(xinm(:,2,xmc))),'-.','linewidth',1.5)
legend({'$\theta - \hat{\theta}$ Ayaw MMAE', '$\theta - \hat{\theta}$ Byaw MMAE' },'interpreter','latex','FontSize',14)
xlabel('Time [s]','FontSize',14)
ylabel('$|\theta-\hat{\theta}|$','interpreter','latex','FontSize',14)
ylim([-beulLim beulLim])
grid on
%

hLim = 200;
ff(2) = figure(2);clf
ff(2).Position = [horix+movex very-0.75*movey 2*figx 2*figy];
subplot(3,2,1)
plot(clock,phinoise-squeeze(x0outm(:,1,xmc)),'linewidth',1.5)
hold on
plot(clock,phinoise-squeeze(x0inm(:,1,xmc)),'linewidth',1.5)
legend({'$\phi - \hat{\phi}$ Ayaw', '$\phi - \hat{\phi}$ Byaw' },'interpreter','latex','FontSize',10)
% ylim(hLim/75*[-eulLim eulLim])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,2)
plot(clock,thetanoise-squeeze(x0outm(:,2,xmc)),'linewidth',1.5)
hold on
plot(clock,thetanoise-squeeze(x0inm(:,2,xmc)),'linewidth',1.5)
legend 'v true' 'u Ayaw'
% ylim(hLim/75*[-eulLim eulLim])
legend({'$\theta - \hat{\theta}$ Ayaw', '$\theta - \hat{\theta}$ Byaw' },'interpreter','latex','FontSize',10)
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,3)
plot(clock,phinoise-squeeze(xgoutm(:,1,xmc)),'linewidth',1.5)
hold on
plot(clock,phinoise-squeeze(xginm(:,1,xmc)),'linewidth',1.5)
legend({'$\phi - \hat{\phi}$ Ayaw Gyro', '$\phi - \hat{\phi}$ Byaw Gyro' },'interpreter','latex','FontSize',10)
% ylim(hLim/75*[-eulLim eulLim])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,4)
plot(clock,thetanoise-squeeze(xgoutm(:,2,xmc)),'linewidth',1.5)
hold on
plot(clock,thetanoise-squeeze(xginm(:,2,xmc)),'linewidth',1.5)
legend({'$\theta - \hat{\theta}$ Ayaw Gyro', '$\theta - \hat{\theta}$ Byaw Gyro' },'interpreter','latex','FontSize',10)
% ylim(hLim/75*[-eulLim eulLim])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,5)
plot(clock,phinoise-squeeze(xaoutm(:,1,xmc)),'linewidth',1.5)
hold on
plot(clock,phinoise-squeeze(xainm(:,1,xmc)),'linewidth',1.5)
legend({'$\phi - \hat{\phi}$ Ayaw Acc', '$\phi - \hat{\phi}$ Byaw Acc' },'interpreter','latex','FontSize',10)
% ylim(hLim/75*[-eulLim eulLim])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,6)
plot(clock,thetanoise-squeeze(xaoutm(:,2,xmc)),'linewidth',1.5)
hold on
plot(clock,thetanoise-squeeze(xainm(:,2,xmc)),'linewidth',1.5)
legend({'$\theta - \hat{\theta}$ Ayaw Acc', '$\theta - \hat{\theta}$ Byaw Acc' },'interpreter','latex','FontSize',10)
% ylim(hLim/75*[-eulLim eulLim])
xlabel 'Time [s]'
xlim([0 120])

%
ff(3) = figure(3);clf
ff(3).Position = [horix+2*movex very figx figy];
h3 = axes;h3.FontSize = 14;hold on
subplot(2,1,1)
plot(clock,windData(:,1),'linewidth',1.5)
hold on
plot(clock,squeeze(x0outm(:,7,xmc)),'linewidth',1.5)
plot(clock,squeeze(xgoutm(:,9,xmc)),'linewidth',1.5)
plot(clock,squeeze(xaoutm(:,9,xmc)),'linewidth',1.5)
legend 'u_w true' 'u_w 0' 'u_w g' 'u_w a'
xlabel 'Time [s]'
ylabel 'u_w [m/s]'
ylim([-windxLim windxLim])

subplot(2,1,2)
plot(clock,windData(:,2),'linewidth',1.5)
hold on
plot(clock,squeeze(x0outm(:,8,xmc)),'linewidth',1.5)
plot(clock,squeeze(xgoutm(:,10,xmc)),'linewidth',1.5)
plot(clock,squeeze(xaoutm(:,10,xmc)),'linewidth',1.5)
legend 'v_w true' 'v_w 0' 'v_w g' 'v_w a'
xlabel 'Time [s]'
ylabel 'v_w [m/s]'
ylim([-windyLim windyLim])


%
ff(4) = figure(4);clf
ff(4).Position = [horix+3*movex very figx figy];
h4 = axes;h4.FontSize = 14;hold on
subplot(2,1,1)
plot(clock,bias_p.*gain_bp,'linewidth',1.5)
hold on
plot(clock,squeeze(xgoutm(:,7,xmc)),'linewidth',1.5)
plot(clock,squeeze(xginm(:,7,xmc)),'-.','linewidth',1.5)
legend 'bias p true' 'bias p out' 'bias p in'
xlabel 'Time [s]'
ylabel 'bias p [rad/s]'
ylim([-beulLim beulLim])

subplot(2,1,2)
plot(clock,bias_q.*gain_bq,'linewidth',1.5)
hold on
plot(clock,squeeze(xgoutm(:,8,xmc)),'linewidth',1.5)
plot(clock,squeeze(xginm(:,8,xmc)),'-.','linewidth',1.5)
legend 'bias q true' 'bias q aut' 'bias q in'
xlabel 'Time [s]'
ylabel 'bias q [rad/s]'
ylim([-beulLim beulLim])

ff(5) = figure(5);clf
ff(5).Position = [horix very-movey figx figy];
h4 = axes;h4.FontSize = 14;hold on
subplot(4,1,1)
plot(clock,bias_p.*gain_bp,'linewidth',1.5)
hold on
plot(clock,squeeze(xgoutm(:,7,xmc)),'linewidth',1.5)
plot(clock,squeeze(xginm(:,7,xmc)),'--','linewidth',1.5,'color',[0.3010 0.7450 0.9330])
legend 'bias p' 'bias p Ayaw Gyro' 'bias p Byaw Gyro'
xlabel 'Time [s]'
ylabel 'bias p [rad/s]'
ylim([-beulLim beulLim])

subplot(4,1,2)
plot(clock,bias_q.*gain_bq,'linewidth',1.5)
hold on
plot(clock,squeeze(xgoutm(:,8,xmc)),'linewidth',1.5)
plot(clock,squeeze(xginm(:,8,xmc)),'-.','linewidth',1.5,'color',[0.3010 0.7450 0.9330])
legend 'bias q' 'bias q Ayaw Gyro' 'bias q Byaw Gyro'
xlabel 'Time [s]'
ylabel 'bias q [rad/s]'
ylim([-beulLim beulLim])

subplot(4,1,3)
plot(clock,bias_ax.*gain_bx,'linewidth',1.5)
hold on
plot(clock,squeeze(xaoutm(:,7,xmc)),'linewidth',1.5)
plot(clock,squeeze(xainm(:,7,xmc)),'-.','linewidth',1.5,'color',[0.3010 0.7450 0.9330])
legend 'bias x' 'bias x Ayaw Acc' 'bias x Byaw Acc'
xlabel 'Time [s]'
ylabel 'bias x [m/s]'
ylim([-beulLim beulLim])

subplot(4,1,4)
plot(clock,bias_ay.*gain_by,'linewidth',1.5)
hold on
plot(clock,squeeze(xaoutm(:,8,xmc)),'linewidth',1.5)
plot(clock,squeeze(xainm(:,8,xmc)),'-.','linewidth',1.5,'color',[0.3010 0.7450 0.9330])
legend 'bias y' 'bias y Ayaw Acc' 'bias y Byaw Acc'
xlabel 'Time [s]'
ylabel 'bias y [m/s]'
ylim([-beulLim beulLim])

%
ff(6) = figure(6);clf
ff(6).Position = [horix+movex very-movey figx figy];
h4 = axes;h4.FontSize = 14;hold on
subplot(2,1,1)
plot(clock(1:200:end),squeeze(P_0outm(1:200:end,1,xmc)),'linewidth',2)
hold on
plot(clock(1:200:end),squeeze(P_goutm(1:200:end,1,xmc)),'linewidth',2)
plot(clock(1:200:end),squeeze(P_aoutm(1:200:end,1,xmc)),'-.','linewidth',2)
plot(clock(1:200:end),squeeze(P_0outm(1:200:end,1,xmc))+squeeze(P_goutm(1:200:end,1,xmc))+squeeze(P_aoutm(1:200:end,1,xmc)),'-.','linewidth',2)
legend 'P Ayaw' 'P Ayaw Gyro' 'P Ayaw Acc' '\Sigma P Ayaw'
xlabel 'Time [s]'
ylabel 'Probability'
ylim([-0.1 1.1])

subplot(2,1,2)
plot(clock(1:200:end),squeeze(P_0inm(1:200:end,1,xmc)),'linewidth',2)
hold on
plot(clock(1:200:end),squeeze(P_ginm(1:200:end,1,xmc)),'linewidth',2)
plot(clock(1:200:end),squeeze(P_ainm(1:200:end,1,xmc)),'-.','linewidth',2)
plot(clock(1:200:end),squeeze(P_0inm(1:200:end,1,xmc))+squeeze(P_ginm(1:200:end,1,xmc))+squeeze(P_ainm(1:200:end,1,xmc)),'-.','linewidth',2)
legend 'P Byaw' 'P Byaw Gyro' 'P Byaw Acc' '\Sigma P Byaw'
xlabel 'Time [s]'
ylabel 'Probability'
ylim([-0.1 1.1])

% ff(7) = figure(7);clf
% ff(7).Position = [horix+2*movex very-movey 1.5*figx 1.5*figy];
% h4 = axes;h4.FontSize = 14;hold on
% subplot(4,1,1)
% plot(clock,bias_p.*gain_bp,'linewidth',1.5)
% hold on
% plot(clock,squeeze(bias_gyroout(:,1,xmc)),'linewidth',1.5)
% plot(clock,squeeze(bias_gyroin(:,1,xmc)),'--','linewidth',1.5,'color',[0.3010 0.7450 0.9330])
% legend 'bias p' 'bias p Ayaw Gyro' 'bias p Byaw Gyro'
% xlabel 'Time [s]'
% ylabel 'bias p [rad/s]'
% ylim([-beulLim beulLim])
%
% subplot(4,1,2)
% plot(clock,bias_q.*gain_bq,'linewidth',1.5)
% hold on
% plot(clock,squeeze(bias_gyroout(:,2,xmc)),'linewidth',1.5)
% plot(clock,squeeze(bias_gyroin(:,2,xmc)),'--','linewidth',1.5,'color',[0.3010 0.7450 0.9330])
% legend 'bias q' 'bias q Ayaw Gyro' 'bias q Byaw Gyro'
% xlabel 'Time [s]'
% ylabel 'bias q [rad/s]'
% ylim([-beulLim beulLim])
%
% subplot(4,1,3)
% plot(clock,bias_ax.*gain_bx,'linewidth',1.5)
% hold on
% plot(clock,squeeze(bias_accout(:,1,xmc)),'linewidth',1.5)
% plot(clock,squeeze(bias_accin(:,1,xmc)),'--','linewidth',1.5,'color',[0.3010 0.7450 0.9330])
% legend 'bias x' 'bias x Ayaw Acc' 'bias x Byaw Acc'
% xlabel 'Time [s]'
% ylabel 'bias x [m/s]'
% ylim([-beulLim beulLim])
%
% subplot(4,1,4)
% plot(clock,bias_ay.*gain_by,'linewidth',1.5)
% hold on
% plot(clock,squeeze(bias_accout(:,2,xmc)),'linewidth',1.5)
% plot(clock,squeeze(bias_accin(:,2,xmc)),'--','linewidth',1.5,'color',[0.3010 0.7450 0.9330])
% legend 'bias y' 'bias y Ayaw Acc' 'bias y Byaw Acc'
% xlabel 'Time [s]'
% ylabel 'bias y [m/s]'
% ylim([-beulLim beulLim])





%
ff(7) = figure(7);clf
ff(7).Position = [horix+2*movex very-movey 2.5*figx 2.5*figy];
h4 = axes;h4.FontSize = 14;hold on
subplot(3,2,1)
hold on
for ii = 1:Nmc
    plot(clock,phinoise-squeeze(x0outm(:,1,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\phi - \hat{\phi}$ Ayaw', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,2)
hold on
for ii = 1:Nmc
    plot(clock,phinoise-squeeze(x0inm(:,1,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\phi - \hat{\phi}$ Byaw', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,3)
hold on
for ii = 1:Nmc
    plot(clock,phinoise-squeeze(xgoutm(:,1,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\phi - \hat{\phi}$ Ayaw Gyro', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,4)
hold on
for ii = 1:Nmc
    plot(clock,phinoise-squeeze(xginm(:,1,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\phi - \hat{\phi}$ Byaw Gyro', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,5)
hold on
for ii = 1:Nmc
    plot(clock,phinoise-squeeze(xaoutm(:,1,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\phi - \hat{\phi}$ Ayaw Acc', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])


subplot(3,2,6)
hold on
for ii = 1:Nmc
    plot(clock,phinoise-squeeze(xainm(:,1,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\phi - \hat{\phi}$ Byaw Acc', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

%
ff(8) = figure(8);clf
ff(8).Position = [horix+3*movex very-movey 2.5*figx 2.5*figy];
h4 = axes;h4.FontSize = 14;hold on
subplot(3,2,1)
hold on
for ii = 1:Nmc
    plot(clock,thetanoise-squeeze(x0outm(:,2,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\theta - \hat{\theta}$ Ayaw', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,2)
hold on
for ii = 1:Nmc
    plot(clock,thetanoise-squeeze(x0inm(:,2,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\theta - \hat{\theta}$ Byaw', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,3)
hold on
for ii = 1:Nmc
    plot(clock,thetanoise-squeeze(xgoutm(:,2,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\theta - \hat{\theta}$ Ayaw Gyro', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,4)
hold on
for ii = 1:Nmc
    plot(clock,thetanoise-squeeze(xginm(:,2,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\theta - \hat{\theta}$ Byaw Gyro', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])

subplot(3,2,5)
hold on
for ii = 1:Nmc
    plot(clock,thetanoise-squeeze(xaoutm(:,2,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\theta - \hat{\theta}$ Ayaw Acc', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])


subplot(3,2,6)
hold on
for ii = 1:Nmc
    plot(clock,thetanoise-squeeze(xainm(:,2,ii)),'linewidth',1.5)
end
legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
    '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
title({'$\theta - \hat{\theta}$ Byaw Acc', },'interpreter','latex','FontSize',14)
% ylim([-0.2 0.2])
xlabel 'Time [s]'
xlim([0 120])


% subplot(2,2,1)
% hold on
% for ii = 1:3:Nmc
% plot(clock,phinoise-squeeze(xgoutm(:,1,ii)),'linewidth',1.5)
% end
% legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
%     '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
% title({'$\phi - \hat{\phi}$ Ayaw Acc', },'interpreter','latex','FontSize',14)
% % ylim([-0.2 0.2])
% xlabel 'Time [s]'
% xlim([0 120])
%
% subplot(2,2,2)
% hold on
% for ii = 1:3:Nmc
% plot(clock,phinoise-squeeze(xginm(:,1,ii)),'linewidth',1.5)
% end
% legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
%     '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
% title({'$\phi - \hat{\phi}$ Byaw Acc'  },'interpreter','latex','FontSize',14)
% xlabel 'Time [s]'
% xlim([0 120])
% % ylim([-0.2 0.2])
%
% subplot(2,2,3)
% hold on
% for ii = 1:3:Nmc
% plot(clock,thetanoise-squeeze(xgoutm(:,2,ii)),'linewidth',1.5)
% end
% legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
%     '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
% title({'$\theta - \hat{\theta}$ Ayaw Acc' },'interpreter','latex','FontSize',14)
% % ylim(hLim/75*[-eulLim eulLim])
% xlabel 'Time [s]'
% xlim([0 120])
% % ylim([-0.2 0.2])
%
%
%
% subplot(2,2,4)
% hold on
% for ii = 1:3:Nmc
% plot(clock,thetanoise-squeeze(xginm(:,2,ii)),'linewidth',1.5)
% end
% legend({'$R_r = 10^{-1},\,R_\psi = 5\times 10^{-2}$','$R_r = 10^{-4},\,R_\psi = 5\times 10^{-5}$',...
%     '$R_r = 10^{-7},\,R_\psi = 5\times 10^{-8}$','$R_r = 10^{-10},\,R_\psi = 5\times 10^{-11}$'},'interpreter','latex','FontSize',10)
% title({'$\theta - \hat{\theta}$ Byaw Acc' },'interpreter','latex','FontSize',14)
% % ylim(hLim/75*[-eulLim eulLim])
% xlabel 'Time [s]'
% xlim([0 120])
% % ylim([-0.2 0.2])





%% Using the integration to estimate the gyro bias

% Integration of pnoise
clear ind_pm ind_qm
max_phi = 0.5;
max_theta = 0.5;
ind_p = 1;
ind_q = 1;
std_p = 0.00;kk_p = 1;kk_q = 1;
bias_pp = 1.0*bias_p;
bias_qq = 1.0*bias_q;
% bias_pp = [zeros(1,2000) 2/22001*(1:22001)]';
thres = 0.01;
xt = [0;0];
T_w = 5;
N_w = T_w/param.dt;

for ii = 1:N
    Fii = [1 psidot(ii)*param.dt;-psidot(ii)*param.dt 1];
    xt = Fii*xt + param.dt*[pnoise(ii)+std_p*randn + bias_pp(ii);...
        qnoise(ii)+std_p*randn + bias_qq(ii)];
    phiest(ii,:) = sum(pnoise(1:ii)+std_p*randn + bias_pp(1:ii))*param.dt;
    thetaest(ii,:) = sum(qnoise(1:ii)+std_p*randn + bias_qq(1:ii))*param.dt;
    bias_pint(ii,:) = sum(pnoise(ind_p:ii)+std_p*randn + bias_pp(ind_p:ii))/(ii-ind_p+1);
    bias_qint(ii,:) = sum(qnoise(ind_q:ii)+std_p*randn + bias_qq(ind_q:ii))/(ii-ind_q+1);
    
    
    
    if abs(phiest(ii,:)-max_phi)<thres||abs(phiest(ii,:)+max_phi)<thres
        ind_pm(kk_p) = ii;
        if kk_p>1
            ind_p = ind_pm(end);
        else
            ind_p = ind_pm(kk_p);
        end
        kk_p = kk_p+1;
    end
    if abs(thetaest(ii,:)-max_theta)<thres*5||abs(thetaest(ii,:)+max_theta)<thres*5
        ind_qm(kk_q) = ii;
        if kk_q>1
            ind_q = ind_qm(end);
        else
            ind_q = ind_qm(kk_q);
        end
        kk_q = kk_q+1;
    end
    
    phiest2(ii,:) = sum(pnoise(1:ii)+std_p*randn + bias_pp(1:ii)-bias_pint(1:ii))*param.dt;
    thetaest2(ii,:) = sum(qnoise(1:ii)+std_p*randn + bias_qq(1:ii)-bias_qint(1:ii))*param.dt;
    
    % Using sliding window
    if ii<N_w
        bias_pw(ii,:) = sum([zeros(N_w-ii,1); pnoise(1:ii)+bias_pp(1:ii)])*param.dt/T_w;
        bias_qw(ii,:) = sum([zeros(N_w-ii,1); qnoise(1:ii)+bias_qq(1:ii)])*param.dt/T_w;
    else
        bias_pw(ii,:) = sum(pnoise(ii-N_w+1:ii)+bias_pp(ii-N_w+1:ii))*param.dt/T_w;
        bias_qw(ii,:) = sum(qnoise(ii-N_w+1:ii)+bias_qq(ii-N_w+1:ii))*param.dt/T_w;
    end
    phiw(ii,:) = sum(pnoise(1:ii)+std_p*randn + bias_pp(1:ii)-bias_pw(1:ii))*param.dt;
    thetaw(ii,:) = sum(qnoise(1:ii)+std_p*randn + bias_qq(1:ii)-bias_qw(1:ii))*param.dt;
end

disp('Done estimate the sensor gyro biases')
%
figure(1);clf
f1 = axes;
f1.FontSize = 18;hold on
subplot(2,1,1)
plot(clock,bias_pp,'-','linewidth',2.5,'color','b')
hold on
plot(clock,bias_pint,'-','linewidth',1.5)
plot(clock,bias_pw,'-.','linewidth',1.5)
xlabel('Time [s]','FontSize',14)
ylabel('$\beta_\phi$','Interpreter','latex','FontSize',18)
legend({'$\beta_\phi$', '$\hat{\beta}_\phi$ int', '$\hat{\beta}_\phi$ fixed-window'},'Interpreter','latex','FontSize',14)
ylim([-0.75 0.75]);
xlim([0 120])

subplot(2,1,2)
plot(clock,bias_qq,'-','linewidth',2.5,'color','b')
hold on
plot(clock,bias_qint,'-','linewidth',1.5)
plot(clock,bias_qw,'-.','linewidth',1.5)
xlabel('Time [s]','FontSize',14)
ylabel('$\beta_\theta$','Interpreter','latex','FontSize',18)
legend({'$\beta_\phi$', '$\hat{\beta}_\phi$ int', '$\hat{\beta}_\theta$ fixed-window'},'Interpreter','latex','FontSize',14)
ylim([-0.75 0.75]);
xlim([0 120])
%
Hfilt = tf(1,[1 1]);
[bf,af] = tfdata(c2d(Hfilt,param.dt,'zoh'),'v');
bias_qfest = filter(bf,af,bias_qint);
bias_qfw = filter(bf,af,bias_qw);
% plot(clock,bias_qfest)

Hfilt = tf(1,[1 1]);
[bf,af] = tfdata(c2d(Hfilt,param.dt,'zoh'),'v');
bias_pfest = filter(bf,af,bias_pint);
% plot(clock,bias_pfest)
bias_pfw = filter(bf,af,bias_pw);

%
for ii = 1:N
    phiest3(ii,:) = sum(pnoise(1:ii)+std_p*randn + bias_pp(1:ii)-bias_pfest(1:ii))*param.dt;
    thetaest3(ii,:) = sum(qnoise(1:ii)+std_p*randn + bias_qq(1:ii)-bias_qfest(1:ii))*param.dt;
end
%
figure(2);clf
plot(clock,bias_pp,'--','linewidth',1.5)
hold on
plot(clock,bias_qq,'--','linewidth',1.5)
plot(clock,bias_qfest)
plot(clock,bias_pfest)
plot(clock,bias_pfw)
plot(clock,bias_qfw)
legend 'p0' 'q0' 'qf' 'pf' 'pfw' 'qfw'

figure(3);clf
plot(phinoise)
hold on
plot(phiest)
plot(phiest3)
plot(phiw)
plot([1 N],[max_phi max_phi])
plot([1 N],[-max_phi -max_phi])

figure(4);clf
plot(thetanoise)
hold on
plot(thetaest)
plot(thetaest3)
plot(thetaw)
plot([1 N],[max_theta max_theta])
plot([1 N],[-max_theta -max_theta])

%% Estimate the acc bias using estimated gyro bias compensation
tic
Qeul = 1e-5*[1 1]; % 1e-7 make no detection, 1e-9 false alarm when acc bias,
Qv = 1e-5*[1 1]; % Should trust more on the euler and less in the bias, no matter the model 1e-5 1e-1
Qx = 1e-5*[1 1];
Qbiasg = 1e-5*[1 1]; % 5 10-5 5 10-5 10-4
Qbiasa = 5*1e-1*[1 1];
Qwind = 1e0*[1 1];
index = index;
close all
clear x7m y7m x8m y8m x9m y9m
% F7 = @(psit, r) eye(10)+param.dt*[  0 r, 0 0, 0 0, 0 0, 0 0;...
%         -r 0, 0 0, 0 0, 0 0, 0 0;...
%         0 param.g(3), -param.lambda1/param.m r, 0 0, 0 0, -param.lambda1/param.m*[cos(psit) sin(psit)] ;...
%        -param.g(3) 0, -r -param.lambda1/param.m, 0 0, 0 0, -param.lambda1/param.m*[-sin(psit) cos(psit)];...
%         0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
%         0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
%         0 0, 0 0, 0 0, 0 0, 0 0;...
%         0 0, 0 0, 0 0, 0 0, 0 0;...
%         0 0, 0 0, 0 0, 0 0, 0 0;...
%         0 0, 0 0, 0 0, 0 0, 0 0];
%
% % x = [phi theta u v x y biasx biasy]
% F8 = @(psit, r) eye(8)+param.dt*[  0 r, 0 0, 0 0, 0 0;...
%     -r 0, 0 0, 0 0, 0 0;...
%     0 param.g(3), 0 r, 0 0, 1 0;...
%    -param.g(3) 0, -r 0, 0 0, 0 1;...
%     0 0, cos(psit) -sin(psit), 0 0, 0 0;...
%     0 0, sin(psit) cos(psit), 0 0, 0 0;...
%     0 0, 0 0, 0 0, 0 0;...
%     0 0, 0 0, 0 0, 0 0];
% % x9 = [phi theta ax ay u v x y biasx biasy]
% F9 = @(psit, r) eye(10) + param.dt*[  0 r, 0 0, 0 0, 0 0, 0 0;...
%     -r 0, 0 0, 0 0, 0 0, 0 0;...
%     0 0, 0 0, 0 0, 0 0, 0 0;...
%     0 0, 0 0, 0 0, 0 0, 0 0;...
%     0 0, 1 0, 0 0, 0 0, 0 0;...
%     0 0, 0 1, 0 0, 0 0, 0 0;...
%     0 0, 0 0, cos(psit) -sin(psit), 0 0, 0 0;...
%     0 0, 0 0, sin(psit) cos(psit), 0 0, 0 0;...
%     0 0, 0 0, 0 0, 0 0, 0 0;...
%     0 0, 0 0, 0 0, 0 0, 0 0];
%
% H71 = @(psit) [0 0, 0 0, 1 0, 0 0, 0 0;...
%        0 0, 0 0, 0 1, 0 0, 0 0;...
%        0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
%        0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0];
% H72 = @(psit) [0 0, param.lambda1/param.m 0, 0 0, 1 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
%        0 0, 0 param.lambda1/param.m, 0 0, 0 1, param.lambda1/param.m*[-sin(psit) cos(psit)]];
% H7 = @(psit) [0 0, 0 0, 1 0, 0 0, 0 0;...
%        0 0, 0 0, 0 1, 0 0, 0 0;...
%        0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
%        0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
%        0 0, param.lambda1/param.m 0, 0 0, 1 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
%        0 0, 0 param.lambda1/param.m, 0 0, 0 1, param.lambda1/param.m*[-sin(psit) cos(psit)]];
% H8 = @(psit) [0 0, 0 0, 1 0, 0 0;...
%     0 0, 0 0, 0 1, 0 0;...
%     0 0, cos(psit) -sin(psit), 0 0, 0 0;...
%     0 0, sin(psit) cos(psit), 0 0, 0 0];
% H91 = @(psit,r) [0 0, 0 0, 0 0, 1 0, 0 0;...
%        0 0, 0 0, 0 0, 0 1, 0 0;...
%        0 0, 0 0, cos(psit) -sin(psit), 0 0, 0 0;...
%        0 0, 0 0, sin(psit) cos(psit), 0 0, 0 0];
% H92 = @(psit,r) [0 param.g(3), -1 0, 0 r, 0 0, 1 0;...
%        -param.g(3) 0, 0 -1, -r 0, 0 0, 0 1];

B7 = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0]';
B8 = param.dt*[1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, -1 0, 0 0, 0 0;...
    0 0, 0 -1, 0 0, 0 0]';
B9 = B7;

G7 = [1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, 0 0, 1 0 0 0;...
    0 0, 0 0, 0 0, 0 1 0 0;...
    0 0, 0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 0, 0 1];
G8 = [1 0, 0 0,0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 1];
G9 = [1 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0;...
    0 0, 0 param.dt, 0 0, 0 0;...
    0 0, param.dt^2/2 0, 0 0, 0 0;...
    0 0, 0 param.dt^2/2, 0 0, 0 0;...
    0 0, 0 0, param.dt 0, 0 0;...
    0 0, 0 0, 0 param.dt, 0 0;...
    0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 1];
G9 = [1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 0, 0 1, 0 0;...
    0 0, 0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 0, 0 1];

Rpos = 9*0.01; Rspeed = 0.01; Reul = 0.01; Racc  = 1e-6;

R71 = diag([Rpos Rpos Rspeed Rspeed]);...
    R72 = diag([Racc Racc]);
R7 = blkdiag(R71,R72);
R8 = R71;
R91 = diag([Rpos Rpos Rspeed Rspeed]);
R92 = 1e2*diag([Racc Racc]);
R9 = blkdiag(R91,R92);

% We use monte carlo simulation
Nmc = 100;
x7m = zeros(N,10,Nmc);
y7m = zeros(N,6,Nmc);
x8m = zeros(N,8,Nmc);
y8m = zeros(N,4,Nmc);
x9m = zeros(N,10,Nmc);
y9m = zeros(N,6,Nmc);
x10m = zeros(N,8,Nmc);
y10m = zeros(N,4,Nmc);

for mc = 1:Nmc
    
    
    bias_poff = 1;bias_qoff = 1;
    
    stateData = stateNon1w1psi;
    
    utrue = stateData(:,13);
    vtrue = stateData(:,14);
    xtrue = stateData(:,1);
    ytrue = stateData(:,2);
    udot = stateData(:,19);
    vdot = stateData(:,20);
    xnoise = xtrue + sqrt(Rpos)*randn(N,1);
    ynoise = ytrue + sqrt(Rpos)*randn(N,1);
    xdotnoise = stateData(:,7) + sqrt(Rspeed)*randn(N,1);
    ydotnoise = stateData(:,8) + sqrt(Rspeed)*randn(N,1);
    
    
    
    %
    Rpsi = 5*1e-5;Rr =1e-4;
    psinoise = stateData(1:end,6) + sqrt(Rpsi)*randn(N,1);
    psidot = stateData(1:end,18) + sqrt(Rr)*randn(N,1); % Similar performance to state 18, r
    %
    
    Q7 = 1e0*1e-3*diag([1e4*Qeul, 1e7*param.dt*Qv, param.dt^2*1e7*Qx, 5*1e-1*Qbiasa, 1e3*Qwind]); % 1e-4 1e-4 1e-8
    Q8 = 1e0*1e-3*diag([1e4*Qeul, 1e7*param.dt*Qv, param.dt^2*1e7*Qx, 5*1e-1*Qbiasa]);
    Q9 = 1e0*1e-6*diag([1e7*Qeul, 1e10*Qv, param.dt*1e8*Qv, param.dt^2*1e8*Qx, 1e1*Qbiasa]);
    
    % Tuning parameter
    Q7 = 1e-3*diag([5*1e4*Qeul, 1e0*Qv, param.dt*1e0*Qx, 5*1e-1*Qbiasa, 1e3*Qwind]); % 1e-4 1e-4 1e-8
    Q8 = 1e-3*diag([5*1e4*Qeul, 1e0*Qv, param.dt*1e0*Qx, 5*1e-1*Qbiasa]);
    Q9 = 1e-6*diag([1e3*Qeul, 1e10*Qv, 1e6*Qv, 1e4*Qx, 1e2*Qbiasa]);
    
    % % Tuning parameter
    % Q7 = 1e-3*diag([1e2*Qeul, 1e0*Qv, 1e0*Qx, 1e-2*Qbiasa, 5*Qwind]); % 1e-4 1e-4 1e-8
    % Q8 = 1e-3*diag([1e2*Qeul, 1e2*Qv, 1e0*Qx, 1e-2*Qbiasa]);
    % Q9 = 1e-6*diag([1e3*Qeul, 1e10*Qv, 1e6*Qv, 1e4*Qx, 1e2*Qbiasa]);
    
    x7 = [zeros(8,1);0*windData(1,1:2)']; x8 = zeros(8,1);x9 = zeros(10,1);x10 = x8;
    P7 = 100*eye(10); P8 = 100*eye(8); P9 = 100*eye(10); P10 = P8;
    for ii = 1:length(clock)
        % Data, Using window and integration method
        u7 = [pnoise(ii)+bias_poff*(bias_pp(ii)-bias_pw(ii)) qnoise(ii)+bias_qoff*(bias_qq(ii)-bias_qw(ii))];
        y7 = [xnoise(ii) ynoise(ii),...
            xdotnoise(ii) ydotnoise(ii),...
            axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii) ];
        u8 = [pnoise(ii)+bias_poff*(bias_pp(ii)-bias_pw(ii)) qnoise(ii)+bias_qoff*(bias_qq(ii)-bias_qw(ii)),...
            axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii) ];
        y8 = [xnoise(ii) ynoise(ii),...
            xdotnoise(ii) ydotnoise(ii)];
        
        u9 = [pnoise(ii)+bias_poff*(bias_pp(ii)-bias_pint(ii)) qnoise(ii)+bias_qoff*(bias_qq(ii)-bias_qint(ii))];
        y9 = [xnoise(ii) ynoise(ii),...
            xdotnoise(ii) ydotnoise(ii),...
            axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii) ];
        u10 = [pnoise(ii)+bias_poff*(bias_pp(ii)-bias_pint(ii)) qnoise(ii)+bias_qoff*(bias_qq(ii)-bias_qint(ii)),...
            axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii) ];
        y10 = [xnoise(ii) ynoise(ii),...
            xdotnoise(ii) ydotnoise(ii)];
        
        psit = psinoise(ii);r = psidot(ii);
        F7ii = eye(10)+param.dt*[  0 r, 0 0, 0 0, 0 0, 0 0;...
            -r 0, 0 0, 0 0, 0 0, 0 0;...
            0 param.g(3), -param.lambda1/param.m r, 0 0, 0 0, -param.lambda1/param.m*[cos(psit) sin(psit)] ;...
            -param.g(3) 0, -r -param.lambda1/param.m, 0 0, 0 0, -param.lambda1/param.m*[-sin(psit) cos(psit)];...
            0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
            0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0];
        H71ii = [0 0, 0 0, 1 0, 0 0, 0 0;...
            0 0, 0 0, 0 1, 0 0, 0 0;...
            0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
            0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0];
        H72ii = [0 0, param.lambda1/param.m 0, 0 0, 1 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
            0 0, 0 param.lambda1/param.m, 0 0, 0 1, param.lambda1/param.m*[-sin(psit) cos(psit)]];
        H7ii = [H71ii;H72ii];
        
        F8ii = eye(8)+param.dt*[  0 r, 0 0, 0 0, 0 0;...
            -r 0, 0 0, 0 0, 0 0;...
            0 param.g(3), 0 r, 0 0, 1 0;...
            -param.g(3) 0, -r 0, 0 0, 0 1;...
            0 0, cos(psit) -sin(psit), 0 0, 0 0;...
            0 0, sin(psit) cos(psit), 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0];
        H8ii = [0 0, 0 0, 1 0, 0 0;...
            0 0, 0 0, 0 1, 0 0;...
            0 0, cos(psit) -sin(psit), 0 0, 0 0;...
            0 0, sin(psit) cos(psit), 0 0, 0 0];
        
        F9ii = eye(10) + param.dt*[  0 r, 0 0, 0 0, 0 0, 0 0;...
            -r 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 1 0, 0 0, 0 0, 0 0;...
            0 0, 0 1, 0 0, 0 0, 0 0;...
            0 0, 0 0, cos(psit) -sin(psit), 0 0, 0 0;...
            0 0, 0 0, sin(psit) cos(psit), 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0];
        H91ii = [0 0, 0 0, 0 0, 1 0, 0 0;...
            0 0, 0 0, 0 0, 0 1, 0 0;...
            0 0, 0 0, cos(psit) -sin(psit), 0 0, 0 0;...
            0 0, 0 0, sin(psit) cos(psit), 0 0, 0 0];
        H92ii = [0 param.g(3), -1 0, 0 r, 0 0, 1 0;...
            -param.g(3) 0, 0 -1, -r 0, 0 0, 0 1];
        H9ii = [H91ii;H92ii];
        
        
        % Prediction step
        x7 = F7ii*x7 + B7*u7';
        P7 = F7ii*P7*F7ii' + G7*Q7*G7';
        
        x8 = F8ii*x8 + B8*u8';
        P8 = F8ii*P8*F8ii' + G8*Q8*G8';
        
        x9 = F7ii*x9 + B7*u9';
        P9 = F7ii*P9*F7ii' + G7*Q7*G7';
        
        x10 = F8ii*x10 + B8*u10';
        P10 = F8ii*P10*F8ii' + G8*Q8*G8';
        
        %     rr7(ii) = rank(obsv(F7ii,H7ii));
        %     rr8(ii) = rank(obsv(F8ii,H8ii));
        %     cond7(ii) = cond(obsv(F7ii,H7ii));
        %     cond8(ii) = cond(obsv(F8ii,H8ii));
        
        P = [0.5 0.5 0.01 0.01 1 1 0.9 0.9 1 1];
        
        %     K9 = place(F9ii,H9ii',P')';
        if mod(ii-1,index) == 0
            % Correction step
            S7 = H7ii*P7*H7ii'+R7;
            K7 = P7*H7ii'/S7;
            x7 = x7 + K7*(y7' - H7ii*x7);
            P7 = P7 - K7*H7ii*P7;
            
            S8 = H8ii*P8*H8ii'+R8;
            K8 = P8*H8ii'/S8;
            x8 = x8 + K8*(y8' - H8ii*x8);
            P8 = P8 - K8*H8ii*P8;
            
            S9 = H7ii*P9*H7ii'+R7;
            K9 = P9*H7ii'/S9;
            x9 = x9 + K9*(y9' - H7ii*x9);
            P9 = P9 - K9*H7ii*P9;
            
            S10 = H8ii*P10*H8ii'+R8;
            K10 = P10*H8ii'/S10;
            x10 = x10 + K10*(y10' - H8ii*x10);
            P10 = P10 - K10*H8ii*P10;
        else
            % Correction step
            S7 = H72ii*P7*H72ii'+R72;
            K7 = P7*H72ii'/S7;
            x7 = x7 + K7*(y7(5:end)' - H72ii*x7);
            P7 = P7 - K7*H72ii*P7;
            
            % Correction step
            S9 = H72ii*P9*H72ii'+R72;
            K9 = P9*H72ii'/S9;
            x9 = x9 + K9*(y9(5:end)' - H72ii*x9);
            P9 = P9 - K9*H72ii*P9;
            
        end
        
        
        x7m(ii,:,mc) = x7;
        y7m(ii,:,mc) = H7ii*x7;
        
        x8m(ii,:,mc) = x8;
        y8m(ii,:,mc) = H8ii*x8;
        
        x9m(ii,:,mc) = x9;
        y9m(ii,:,mc) = H9ii*x9;
        
        x10m(ii,:,mc) = x10;
        y10m(ii,:,mc) = H8ii*x10;
    end
    disp('Done compute acc bias ')
    
    toc
    tic
    mc
end

% PLOTING THE RESULT OBTAIN FROM THE RUNS

baccLim = 1;windxLim = 2;windyLim = 2;uLim = 3;vLim = 3;
eulLim = 1;beulLim = 1;vxLim = 2;vyLim = 2;windxLim = 4;windyLim = 2;
horix = -1900 ;very = 700;movex = 450; movey = 450; figx = 420; figy = 320;%560 420

if Nmc == 1
    xmc = 1;
else
    xmc = randi(Nmc,1,1);
end

ff(1) = figure(1);clf
ff(1).Position = [horix very figx figy];
subplot(2,1,1)
plot(clock,phinoise,'linewidth',1.5)
hold on
plot(clock,squeeze(x7m(:,1,xmc)),'-','linewidth',1.5)
plot(clock,squeeze(x8m(:,1,xmc)),'-.','linewidth',1.5,'color','g')
plot(clock,squeeze(x9m(:,1,xmc)),'-','linewidth',1.5)

ylim([-eulLim eulLim])

legend 'phi true' 'phi Ayaw Acc' 'phi Byaw Acc' 'theta Cyaw'
subplot(2,1,2)
plot(clock,thetanoise,'linewidth',1.5)
hold on
plot(clock,squeeze(x7m(:,2,xmc)),'-','linewidth',1.5)
plot(clock,squeeze(x8m(:,2,xmc)),'-.','linewidth',1.5,'color','g')
plot(clock,squeeze(x9m(:,2,xmc)),'-','linewidth',1.5)
legend 'theta true' 'theta Ayaw Acc' 'theta Byaw Acc' 'theta Cyaw'
ylim([-eulLim eulLim])

ff(2) = figure(2);clf
ff(2).Position = [horix+movex very figx figy];
subplot(2,1,1)
plot(clock,utrue,'linewidth',1.5)
hold on
plot(clock,squeeze(x7m(:,3,xmc)),'-','linewidth',1.5)
plot(clock,squeeze(x8m(:,3,xmc)),'-.','linewidth',1.5,'color','g')
plot(clock,squeeze(x9m(:,5,xmc)),'-','linewidth',1.5)
legend 'u true' 'u Ayaw Acc' 'u Byaw Acc'
ylim([-uLim uLim])

subplot(2,1,2)
plot(clock,vtrue,'linewidth',1.5)
hold on
plot(clock,squeeze(x7m(:,4,xmc)),'-','linewidth',1.5)
plot(clock,squeeze(x8m(:,4,xmc)),'-.','linewidth',1.5,'color','g')
plot(clock,squeeze(x9m(:,6,xmc)),'-','linewidth',1.5)
legend 'v true'  'v Ayaw Acc' 'v Byaw Acc'
ylim([-vLim vLim])

%
ff(3) = figure(3);clf
ff(3).Position = [horix+2*movex very figx figy];
subplot(2,1,1)
plot(clock,xtrue,'linewidth',1.5)
hold on
plot(clock,squeeze(x7m(:,5,xmc)),'-','linewidth',1.5)
plot(clock,squeeze(x8m(:,5,xmc)),'-.','linewidth',1.5,'color','g')
plot(clock,squeeze(x9m(:,7,xmc)),'-','linewidth',1.5)
legend 'x true' 'x Ayaw Acc' 'x Byaw Acc'
% ylim([-uLim uLim])

subplot(2,1,2)
plot(clock,ytrue,'linewidth',1.5)
hold on
plot(clock,squeeze(x7m(:,6,xmc)),'-','linewidth',1.5)
plot(clock,squeeze(x8m(:,6,xmc)),'-.','linewidth',1.5,'color','g')
plot(clock,squeeze(x9m(:,8,xmc)),'-','linewidth',1.5)
legend 'y true'  'y Ayaw Acc' 'y Byaw Acc'
% ylim([-vLim vLim])


bias_x7 = squeeze(x7m(:,7,:));
bias_y7 = squeeze(x7m(:,8,:));
bias_x8 = squeeze(x8m(:,7,:));
bias_y8 = squeeze(x8m(:,8,:));
m_bias_x7 = mean(bias_x7,2);
m_bias_y7 = mean(bias_y7,2);
m_bias_x8 = mean(bias_x8,2);
m_bias_y8 = mean(bias_y8,2);

std_bias_x7 = std(bias_x7,0,2);
std_bias_y7 = std(bias_y7,0,2);
std_bias_x8 = std(bias_x8,0,2);
std_bias_y8 = std(bias_y8,0,2);

bias_x9 = squeeze(x9m(:,7,:));
bias_y9 = squeeze(x9m(:,8,:));
bias_x10 = squeeze(x10m(:,7,:));
bias_y10 = squeeze(x10m(:,8,:));
m_bias_x9 = mean(bias_x9,2);
m_bias_y9 = mean(bias_y9,2);
m_bias_x10 = mean(bias_x10,2);
m_bias_y10 = mean(bias_y10,2);

std_bias_x9 = std(bias_x9,0,2);
std_bias_y9 = std(bias_y9,0,2);
std_bias_x10 = std(bias_x10,0,2);
std_bias_y10 = std(bias_y10,0,2);

ff(4) = figure(4);clf
ff(4).Position = [horix+3*movex very figx figy];
subplot(2,1,1)
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,m_bias_x7,'-','linewidth',1.5,'color','r')
plot(clock,m_bias_x8,'-.','linewidth',1.5,'color','g')
plot(clock,m_bias_x7+std_bias_x7,'-','linewidth',1,'color','r')
plot(clock,m_bias_x7-std_bias_x7,'-','linewidth',1,'color','r')
plot(clock,m_bias_x8+std_bias_x8,'--','linewidth',1,'color','g')
plot(clock,m_bias_x8-std_bias_x8,'--','linewidth',1,'color','g')
% plot(clock,x9m(:,9),'-','linewidth',1.5)
legend({'$\beta_x$', '$\hat{\beta}_x$ Ayaw Acc', '$\hat{\beta}_x$ Byaw Acc'},'Interpreter','latex','Location','ne')
grid on
xlabel 'Time [s]'
% legend 'bias ax true' 'bias ax aout' 'bias ax aout+non'
ylim(0.75*[-baccLim baccLim]+mean(bias_ax))
title 'Accelerometer bias estimates with moving window gyro bias estimate'

subplot(2,1,2)
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,m_bias_y7,'-','linewidth',1.5,'color','r')
plot(clock,m_bias_y8,'-.','linewidth',1.5,'color','g')
plot(clock,m_bias_y7+std_bias_y7,'-','linewidth',1,'color','r')
plot(clock,m_bias_y7-std_bias_y7,'-','linewidth',1,'color','r')
plot(clock,m_bias_y8+std_bias_y8,'--','linewidth',1,'color','g')
plot(clock,m_bias_y8-std_bias_y8,'--','linewidth',1,'color','g')
legend({'$\beta_y$', '$\hat{\beta}_y$ Ayaw Acc', '$\hat{\beta}_y$ Byaw Acc'},'Interpreter','latex')
grid on
xlabel 'Time [s]'
ylim(0.75*[-baccLim baccLim]+mean(bias_ay))

ff(5) = figure(5);clf
ff(5).Position = [horix very-movey figx figy];
subplot(2,1,1)
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,m_bias_x9,'-','linewidth',1.5,'color','r')
plot(clock,m_bias_x10,'-.','linewidth',1.5,'color','g')
plot(clock,m_bias_x9+std_bias_x9,'-','linewidth',1,'color','r')
plot(clock,m_bias_x9-std_bias_x9,'-','linewidth',1,'color','r')
plot(clock,m_bias_x10+std_bias_x10,'--','linewidth',1,'color','g')
plot(clock,m_bias_x10-std_bias_x10,'--','linewidth',1,'color','g')

legend({'$\beta_x$', '$\hat{\beta}_x$ Ayaw Acc', '$\hat{\beta}_x$ Byaw Acc'},'Interpreter','latex','Location','ne')
grid on
xlabel 'Time [s]'
ylim(0.75*[-baccLim baccLim]+mean(bias_ax))
title 'Accelerometer bias estimates with integrating gyro bias estimate'

subplot(2,1,2)
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,m_bias_y9,'-','linewidth',1.5,'color','r')
plot(clock,m_bias_y10,'-.','linewidth',1.5,'color','g')
plot(clock,m_bias_y9+std_bias_y9,'-','linewidth',1,'color','r')
plot(clock,m_bias_y9-std_bias_y9,'-','linewidth',1,'color','r')
plot(clock,m_bias_y10+std_bias_y10,'--','linewidth',1,'color','g')
plot(clock,m_bias_y10-std_bias_y10,'--','linewidth',1,'color','g')
legend({'$\beta_y$', '$\hat{\beta}_y$ Ayaw Acc', '$\hat{\beta}_y$ Byaw Acc'},'Interpreter','latex')
grid on
xlabel 'Time [s]'
ylim(0.75*[-baccLim baccLim]+mean(bias_ay))

%
phi_x7 = squeeze(x7m(:,1,:));
theta_x7 = squeeze(x7m(:,2,:));
phi_x8 = squeeze(x8m(:,1,:));
theta_x8 = squeeze(x8m(:,2,:));
m_phi_x7 = mean(phi_x7,2);
m_theta_x7 = mean(theta_x7,2);
m_phi_x8 = mean(phi_x8,2);
m_theta_x8 = mean(theta_x8,2);

std_phi_x7 = std(phi_x7,0,2);
std_theta_x7 = std(theta_x7,0,2);
std_phi_x8 = std(phi_x8,0,2);
std_theta_x8 = std(theta_x8,0,2);

phi_x9 = squeeze(x9m(:,1,:));
theta_x9 = squeeze(x9m(:,2,:));
phi_x10 = squeeze(x10m(:,1,:));
theta_x10 = squeeze(x10m(:,2,:));
m_phi_x9 = mean(phi_x9,2);
m_theta_x9 = mean(theta_x9,2);
m_phi_x10 = mean(phi_x10,2);
m_theta_x10 = mean(theta_x10,2);

std_phi_x9 = std(phi_x9,0,2);
std_theta_x9 = std(theta_x9,0,2);
std_phi_x10 = std(phi_x10,0,2);
std_theta_x10 = std(theta_x10,0,2);

ff(6) = figure(6);clf
ff(6).Position = [horix+movex very-movey 1.5*figx 1.5*figy];
subplot(2,1,1)
hold on
plot(clock,stateData(:,4)-m_phi_x7,'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,4)-m_phi_x8,'-.','linewidth',1.5,'color','g')
% plot(clock,stateData(:,4)-m_phi_x7+std_phi_x7,'--','linewidth',1,'color','r')
% plot(clock,stateData(:,4)-m_phi_x8+std_phi_x8,'--','linewidth',1,'color','g')
% plot(clock,stateData(:,4)-m_phi_x7-std_phi_x7,'--','linewidth',1,'color','r')
% plot(clock,stateData(:,4)-m_phi_x8-std_phi_x8,'--','linewidth',1,'color','g')
legend({'$\phi-\hat{\phi}$ Ayaw Acc','$\phi-\hat{\phi}$ Byaw Acc'},'interpreter','latex','FontSize',14)
xlabel('Time [s]','FontSize',14)
ylim(0.25*[-eulLim eulLim])
grid on

subplot(2,1,2)
plot(clock,stateData(:,5)-m_theta_x7,'-','linewidth',1.5,'color','r')
hold on
plot(clock,stateData(:,5)-m_theta_x8,'-.','linewidth',1,'color','g')

% plot(clock,stateData(:,5)-m_theta_x7+std_theta_x7,'--','linewidth',1,'color','r')
% plot(clock,stateData(:,5)-m_theta_x7-std_theta_x7,'--','linewidth',1,'color','r')
% plot(clock,stateData(:,5)-m_theta_x8+std_theta_x8,'--','linewidth',1,'color','g')
% plot(clock,stateData(:,5)-m_theta_x8-std_theta_x8,'--','linewidth',1,'color','g')
legend({'$\theta-\hat{\theta}$ Ayaw Acc','$\theta-\hat{\theta}$ Byaw Acc'},'interpreter','latex','FontSize',14)
grid on
xlabel('Time [s]','FontSize',14)
ylim(0.25*[-eulLim eulLim])


ff(7) = figure(7);clf
ff(7).Position = [horix+2*movex very-movey figx figy];
subplot(2,1,1)
subplot(2,1,1)
hold on
plot(clock,stateData(:,4)-m_phi_x9,'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,4)-m_phi_x10,'-.','linewidth',1.5,'color','g')
plot(clock,stateData(:,4)-m_phi_x9+std_phi_x9,'--','linewidth',1,'color','r')
plot(clock,stateData(:,4)-m_phi_x10+std_phi_x10,'--','linewidth',1,'color','g')
plot(clock,stateData(:,4)-m_phi_x9-std_phi_x9,'--','linewidth',1,'color','r')
plot(clock,stateData(:,4)-m_phi_x10-std_phi_x10,'--','linewidth',1,'color','g')
legend({'$\phi-\hat{\phi}$ Ayaw Acc','$\phi-\hat{\phi}$ Byaw Acc'},'interpreter','latex','FontSize',14)
xlabel('Time [s]','FontSize',14)
ylim([-eulLim eulLim])
grid on

subplot(2,1,2)
hold on
plot(clock,stateData(:,5)-m_theta_x9,'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,5)-m_theta_x10,'-.','linewidth',1.5,'color','g')
plot(clock,stateData(:,5)-m_theta_x9+std_theta_x9,'-','linewidth',1,'color','r')
plot(clock,stateData(:,5)-m_theta_x10+std_theta_x10,'--','linewidth',1,'color','g')
plot(clock,stateData(:,5)-m_theta_x9-std_theta_x9,'--','linewidth',1,'color','r')
plot(clock,stateData(:,5)-m_theta_x10-std_theta_x10,'--','linewidth',1,'color','g')
legend({'$\theta-\hat{\theta}$ Ayaw Acc','$\theta-\hat{\theta}$ Byaw Acc'},'interpreter','latex','FontSize',14)
grid on
xlabel('Time [s]','FontSize',14)
ylim([-eulLim eulLim])
%

tic
ff(8) = figure(8);clf
title('Sensor bias estimation','FontSize',14)
ff(8).Position = [horix+3*movex very-movey figx figy];
ff(1,1) = subplot(2,2,1);ff(1,1).FontSize = 12;hold on
plot(clock,bias_p,'linewidth',2)
hold on
plot(clock,bias_pw,'-.','linewidth',1.5)
ylabel({'$\beta_\phi$'},'interpreter','latex','FontSize',14)
xlabel({'time [s]'},'interpreter','latex','FontSize',14)
legend({'$\beta_\phi$','$\hat{\beta}_\phi$'},'interpreter','latex','FontSize',12,'Location','se')
ylim(0.5*[-0 eulLim])
xlim([0 120])
grid on

ff(1,2) = subplot(2,2,2);ff(1,2).FontSize = 12;hold on
plot(clock,bias_q,'linewidth',2)
hold on
plot(clock,bias_qw,'-.','linewidth',1.5)
ylabel({'$\beta_\theta$'},'interpreter','latex','FontSize',14)
xlabel({'time [s]'},'interpreter','latex','FontSize',14)
legend({'$\beta_\theta$','$\hat{\beta}_\theta$'},'interpreter','latex','FontSize',12,'Location','ne')
ylim(0.5*[-eulLim 0])
xlim([0 120])
grid on

ff(2,1) = subplot(2,2,3);ff(2,1).FontSize = 12;hold on
plot(clock,bias_ax,'linewidth',2)
hold on
plot(clock,squeeze(x8m(:,7,xmc)),'-.','linewidth',1.5)
ylabel({'$\beta_x$'},'interpreter','latex','FontSize',14)
xlabel({'time [s]'},'interpreter','latex','FontSize',14)
legend({'$\beta_x$','$\hat{\beta}_x$'},'interpreter','latex','FontSize',12,'Location','se')
ylim([-1 1])
xlim([0 120])
grid on

ff(2,2) = subplot(2,2,4);ff(2,2).FontSize = 12;hold on
plot(clock,bias_ay,'linewidth',2)
hold on
plot(clock,squeeze(x8m(:,8,xmc)),'-.','linewidth',1.5)
ylabel({'$\beta_y$'},'interpreter','latex','FontSize',14)
xlabel({'time [s]'},'interpreter','latex','FontSize',14)
legend({'$\beta_y$','$\hat{\beta}_y$'},'interpreter','latex','FontSize',12,'Location','ne')
ylim([-1 1])
xlim([0 120])
grid on


[diag(P7) [diag(P8);0;0] diag(P9)]
toc





%% Using the bias window+full model of acc and gyro bias
Qeul = 1e-5*[1 1]; % 1e-7 make no detection, 1e-9 false alarm when acc bias,
Qv = 1e-5*[1 1]; % Should trust more on the euler and less in the bias, no matter the model 1e-5 1e-1
Qx = 1e-5*[1 1];
Qbiasg = 1e-5*[1 1]; % 5 10-5 5 10-5 10-4
Qbiasa = 5*1e-3*[1 1];
Qwind = 1e0*[1 1];
index = index;

bias_poff = 1;bias_qoff = 1;
reducedNoise = 1e-12;
Rpos = 9*0.01*reducedNoise; Rspeed = 0.01*reducedNoise;
Reul = 1e0; Racc  = 1e-6*reducedNoise;Rgyro = 1e-6*reducedNoise;
close all
clear x7m y7m x8m y8m x9m y9m

% 7: aout+no noise yaw, 8:ain+no noise yaw, 9: aout+noisy yaw, 10 ain+noisy
% yaw
B7 = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0, 0 0]';
B8 = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, -1 0, 0 0, 0 0, 0 0;...
    0 0, 0 -1, 0 0, 0 0, 0 0]';
B9 = B7;
B10 = B8;
G7 = eye(12);G8 = eye(10);G9 = G7;G10 = G8;

R71 = reducedNoise*diag([Rpos Rpos Rspeed Rspeed Reul Reul])+1e-8*eye(6);...
    R72 = reducedNoise*diag([Racc Racc])+1e-8*eye(2) ;
R7 = blkdiag(R71,R72);
R8 = reducedNoise*diag([Rpos Rpos Rspeed Rspeed Reul Reul])+1e-8*eye(6);
R91 = reducedNoise*diag([Rpos Rpos Rspeed Rspeed Reul Reul])+1e-8*eye(6);
R92 = reducedNoise*diag([Racc Racc])+1e-8*eye(2);
R9 = reducedNoise*blkdiag(R91,R92);

% We use monte carlo simulation
Nmc = 100;
x7m = zeros(N,12,Nmc);
y7m = zeros(N,8,Nmc);
x8m = zeros(N,10,Nmc);
y8m = zeros(N,6,Nmc);
x9m = zeros(N,12,Nmc);
y9m = zeros(N,8,Nmc);
x10m = zeros(N,10,Nmc);
y10m = zeros(N,6,Nmc);
bias_pm = zeros(N,1,Nmc);
bias_qm = zeros(N,1,Nmc);
T_w = 1;
N_w = T_w/param.dt;
IMUdata = IMU1w1psi;
stateData = stateNon1w1psi;
windData = wind1w1psi;

utrue = stateData(:,13);
vtrue = stateData(:,14);
xtrue = stateData(:,1);
ytrue = stateData(:,2);
udot = stateData(:,19);
vdot = stateData(:,20);

tic
for mc = 1:Nmc
    
    
    
    xnoise = xtrue + sqrt(Rpos)*randn(N,1);
    ynoise = ytrue + sqrt(Rpos)*randn(N,1);
    xdotnoise = stateData(:,7) + sqrt(Rspeed)*randn(N,1);
    ydotnoise = stateData(:,8) + sqrt(Rspeed)*randn(N,1);
    
    pnoise = IMUdata(:,4) + sqrt(Rgyro)*randn(N,1);
    qnoise = IMUdata(:,5) + sqrt(Rgyro)*randn(N,1);
    rnoise = IMUdata(:,6) + sqrt(Rgyro)*randn(N,1);
    axnoise = -IMUdata(:,1) + sqrt(Racc)*randn(N,1);
    aynoise = -IMUdata(:,2) + sqrt(Racc)*randn(N,1);
    
    Rpsi = 1e-4;Rr =1e-6;
    psifree = stateData(1:end,6);
    psidotfree = stateData(1:end,18);
    psinoise = psifree + sqrt(Rpsi)*randn(N,1);
    psidotnoise = psidotfree + sqrt(Rr)*randn(N,1); % Similar performance to state 18, r
    
    Q7 = reducedNoise*1e0*1e-3*diag([1e4*Qeul, 1e7*param.dt*Qv, param.dt^2*1e7*Qx, 5*1e-1*Qbiasa, Qbiasg, 1e3*Qwind]); % 1e-4 1e-4 1e-8
    Q8 = reducedNoise*1e0*1e-3*diag([1e4*Qeul, 1e10*param.dt*Qv, param.dt^2*1e10*Qx, 5*1e-1*Qbiasa, Qbiasg]);
    Q9 = reducedNoise*1e0*1e-6*diag([1e7*Qeul, 1e10*Qv, param.dt*1e8*Qv, param.dt^2*1e8*Qx, 1e1*Qbiasa]);
    
    % Tuning parameter
    Q7 = reducedNoise*1e-3*diag([5*1e4*Qeul, 1e0*Qv, param.dt*1e0*Qx, 5*1e0*Qbiasa, 1e2*Qbiasg, 1e3*Qwind]); % 1e-4 1e-4 1e-8
    Q8 = reducedNoise*1e-3*diag([5*1e4*Qeul, 1e3*Qv, param.dt*1e5*Qx, 5*1e0*Qbiasa, 1e2*Qbiasg]);
    Q9 = reducedNoise*1e-6*diag([1e3*Qeul, 1e10*Qv, 1e6*Qv, 1e4*Qx, 1e2*Qbiasa, Qbiasg]);
    
    
    % Integration of pnoise
    clear ind_pm ind_qm
    bias_pp = 1.0*bias_p;
    bias_qq = 1.0*bias_q;
    
    x7 = [zeros(10,1);0*windData(1,1:2)']; x8 = zeros(10,1);x9 = zeros(12,1);x10 = x8;
    P7 = 100*eye(12); P8 = 100*eye(10); P9 = 100*eye(12); P10 = P8;
    for ii = 1:length(clock)
        
        % Compute the bias from window method
        % Using sliding window
        if mod(ii-1,index) == 0 % same sampling frequency with GPS
            if ii<N_w
                bias_pwii = sum([zeros(N_w-ii,1); pnoise(1:ii)+bias_pp(1:ii)+sqrt(Rgyro)*randn(ii,1)])*param.dt/T_w;
                bias_qwii = sum([zeros(N_w-ii,1); qnoise(1:ii)+bias_qq(1:ii)+sqrt(Rgyro)*randn(ii,1)])*param.dt/T_w;
            else
                bias_pwii = sum(pnoise(ii-N_w+1:ii)+bias_pp(ii-N_w+1:ii)+sqrt(Rgyro)*randn(N_w,1))*param.dt/T_w;
                bias_qwii = sum(qnoise(ii-N_w+1:ii)+bias_qq(ii-N_w+1:ii)+sqrt(Rgyro)*randn(N_w,1))*param.dt/T_w;
            end
        end
        
        % Data, Using window and integration method
        u7 = [pnoise(ii)+bias_pp(ii) qnoise(ii)+bias_qq(ii)];
        y7 = [xnoise(ii) ynoise(ii),...
            xdotnoise(ii) ydotnoise(ii),...
            bias_pwii bias_qwii,...
            axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii)];
        
        u8 = [pnoise(ii)+bias_pp(ii) qnoise(ii)+bias_qq(ii),...
            axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii) ];
        y8 = [xnoise(ii) ynoise(ii),...
            xdotnoise(ii) ydotnoise(ii),...
            bias_pwii bias_qwii];
        
        psifreet = psifree(ii); rfreet = psidotfree(ii);
        psinoiset = psinoise(ii); rnoiset = psidotnoise(ii);
        
        % F7ii is noise free yaw rate and angle
        F7ii = eye(12)+param.dt*[  0 rfreet, 0 0, 0 0, 0 0, -1 0, 0 0;...
            -rfreet 0, 0 0, 0 0, 0 0, 0 -1, 0 0;...
            0 param.g(3), -param.lambda1/param.m rfreet, 0 0, 0 0, 0 0, -param.lambda1/param.m*[cos(psifreet) sin(psifreet)];...
            -param.g(3) 0, -rfreet -param.lambda1/param.m, 0 0, 0 0, 0 0, -param.lambda1/param.m*[-sin(psifreet) cos(psifreet)];...
            0 0, cos(psifreet) -sin(psifreet), 0 0, 0 0, 0 0, 0 0;...
            0 0, sin(psifreet) cos(psifreet), 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0];
        
        H71ii = [0 0, 0 0, 1 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 1, 0 0, 0 0, 0 0;...
            0 0, cos(psifreet) -sin(psifreet), 0 0, 0 0, 0 0, 0 0;...
            0 0, sin(psifreet) cos(psifreet), 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 1 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 1, 0 0];
        
        H72ii = [0 0, param.lambda1/param.m 0, 0 0, 1 0, 0 0, param.lambda1/param.m*[cos(psifreet) sin(psifreet)];...
            0 0, 0 param.lambda1/param.m, 0 0, 0 1, 0 0, param.lambda1/param.m*[-sin(psifreet) cos(psifreet)]];
        
        H7ii = [H71ii;H72ii];
        
        F8ii = eye(10)+param.dt*[  0 rfreet, 0 0, 0 0, 0 0, -1 0;...
            -rfreet 0, 0 0, 0 0, 0 0, 0 -1;...
            0 param.g(3), 0 rfreet, 0 0, 1 0, 0 0;...
            -param.g(3) 0, -rfreet 0, 0 0, 0 1, 0 0;...
            0 0, cos(psifreet) -sin(psifreet), 0 0, 0 0, 0 0;...
            0 0, sin(psifreet) cos(psifreet), 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0];
        H81ii = [0 0, 0 0, 1 0, 0 0, 0 0;...
            0 0, 0 0, 0 1, 0 0, 0 0;...
            0 0, cos(psifreet) -sin(psifreet), 0 0, 0 0, 0 0;...
            0 0, sin(psifreet) cos(psifreet), 0 0, 0 0, 0 0];
        H82ii = [0 0, 0 0, 0 0, 0 0, 1 0;...
            0 0, 0 0, 0 0, 0 0, 0 1];
        H8ii = [H81ii; H82ii];
        
        F9ii = eye(12)+param.dt*[  0 rnoiset, 0 0, 0 0, 0 0, -1 0, 0 0;...
            -rnoiset 0, 0 0, 0 0, 0 0, 0 -1, 0 0;...
            0 param.g(3), -param.lambda1/param.m rnoiset, 0 0, 0 0, 0 0, -param.lambda1/param.m*[cos(psinoiset) sin(psinoiset)];...
            -param.g(3) 0, -rnoiset -param.lambda1/param.m, 0 0, 0 0, 0 0, -param.lambda1/param.m*[-sin(psinoiset) cos(psinoiset)];...
            0 0, cos(psinoiset) -sin(psinoiset), 0 0, 0 0, 0 0, 0 0;...
            0 0, sin(psinoiset) cos(psinoiset), 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0, 0 0];
        H91ii = [0 0, 0 0, 1 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 1, 0 0, 0 0, 0 0;...
            0 0, cos(psinoiset) -sin(psinoiset), 0 0, 0 0, 0 0, 0 0;...
            0 0, sin(psinoiset) cos(psinoiset), 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 1 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 1, 0 0];
        H92ii = [0 0, param.lambda1/param.m 0, 0 0, 1 0, 0 0, param.lambda1/param.m*[cos(psinoiset) sin(psinoiset)];...
            0 0, 0 param.lambda1/param.m, 0 0, 0 1, 0 0, param.lambda1/param.m*[-sin(psinoiset) cos(psinoiset)]];
        H9ii = [H91ii;H92ii];
        
        F10ii = eye(10)+param.dt*[  0 rnoiset, 0 0, 0 0, 0 0, -1 0;...
            -rnoiset 0, 0 0, 0 0, 0 0, 0 -1;...
            0 param.g(3), 0 rnoiset, 0 0, 1 0, 0 0;...
            -param.g(3) 0, -rnoiset 0, 0 0, 0 1, 0 0;...
            0 0, cos(psinoiset) -sin(psinoiset), 0 0, 0 0, 0 0;...
            0 0, sin(psinoiset) cos(psinoiset), 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0;...
            0 0, 0 0, 0 0, 0 0, 0 0];
        H101ii = [0 0, 0 0, 1 0, 0 0, 0 0;...
            0 0, 0 0, 0 1, 0 0, 0 0;...
            0 0, cos(psinoiset) -sin(psinoiset), 0 0, 0 0, 0 0;...
            0 0, sin(psinoiset) cos(psinoiset), 0 0, 0 0, 0 0];
        H102ii = [0 0, 0 0, 0 0, 0 0, 1 0;...
            0 0, 0 0, 0 0, 0 0, 0 1];
        H10ii = [H101ii; H102ii];
        
        
        % Prediction step
        x7 = F7ii*x7 + B7*u7';
        P7 = F7ii*P7*F7ii' + G7*Q7*G7';

        x8 = F8ii*x8 + B8*u8';
        P8 = F8ii*P8*F8ii' + G8*Q8*G8';

        x9 = F9ii*x9 + B7*u7';
        P9 = F9ii*P9*F9ii' + G7*Q7*G7';

        x10 = F10ii*x10 + B8*u8';
        P10 = F10ii*P10*F10ii' + G8*Q8*G8';
        
        
        if mod(ii-1,index) == 0
            
            % Correction step
            S7 = H7ii*P7*H7ii'+R7;
            K7 = P7*H7ii'/S7;
            x7 = x7 + K7*(y7' - H7ii*x7);
            P7 = P7 - K7*H7ii*P7;
            
            S8 = H8ii*P8*H8ii'+R8;
            K8 = P8*H8ii'/S8;
            x8 = x8 + K8*(y8' - H8ii*x8);
            P8 = P8 - K8*H8ii*P8;
            
            S9 = H9ii*P9*H9ii'+R7;
            K9 = P9*H9ii'/S9;
            x9 = x9 + K9*(y7' - H9ii*x9);
            P9 = P9 - K9*H9ii*P9;
            
            S10 = H10ii*P10*H10ii'+R8;
            K10 = P10*H10ii'/S10;
            x10 = x10 + K10*(y8' - H10ii*x10);
            P10 = P10 - K10*H10ii*P10;
        else
            % Correction step
            S7 = H72ii*P7*H72ii'+R72;
            K7 = P7*H72ii'/S7;
            x7 = x7 + K7*(y7(7:end)' - H72ii*x7);
            P7 = P7 - K7*H72ii*P7;
            
            %         S8 = H82ii*P8*H82ii'+R8(5:end,5:end);
            %         K8 = P8*H82ii'/S8;
            %         x8 = x8 + K8*(y8(5:end)' - H82ii*x8);
            %         P8 = P8 - K8*H82ii*P8;
            
            % Correction step
            S9 = H92ii*P9*H92ii'+R72;
            K9 = P9*H92ii'/S9;
            x9 = x9 + K9*(y7(7:end)' - H92ii*x9);
            P9 = P9 - K9*H92ii*P9;
            
            %         S10 = H102ii*P10*H102ii'+R8(5:end,5:end);
            %         K10 = P10*H102ii'/S10;
            %         x10 = x10 + K10*(y8(5:end)' - H102ii*x10);
            %         P10 = P10 - K10*H102ii*P10;
        end
        
        
        x7m(ii,:,mc) = x7;
        y7m(ii,:,mc) = H7ii*x7;
        
        x8m(ii,:,mc) = x8;
        y8m(ii,:,mc) = H8ii*x8;
        
        x9m(ii,:,mc) = x9;
        y9m(ii,:,mc) = H9ii*x9;
        
        x10m(ii,:,mc) = x10;
        y10m(ii,:,mc) = H10ii*x10;
        
        bias_pm(ii,:,mc) = bias_pwii;
        bias_qm(ii,:,mc) = bias_qwii;
    end
    disp('Done compute acc bias ')
    
    mc
    toc
    tic
    
end

%%

% PLOTING THE RESULT OBTAIN FROM THE RUNS

baccLim = 1;windxLim = 2;windyLim = 2;uLim = 3;vLim = 3;
eulLim = 1;beulLim = 1;vxLim = 2;vyLim = 2;windxLim = 4;windyLim = 2;
horix = -1900 ;very = 700;movex = 450; movey = 450; figx = 420; figy = 320;%560 420

if Nmc == 1
    xmc = 1;
else
    xmc = randi(Nmc,1,1);
end

ff(1) = figure(1);clf
ff(1).Position = [horix very figx figy];
subplot(2,1,1)
plot(clock,phinoise,'linewidth',1.5)
hold on
plot(clock,squeeze(x7m(:,1,xmc)),'-','linewidth',1.5)
plot(clock,squeeze(x8m(:,1,xmc)),'-.','linewidth',1.5,'color','g')
plot(clock,squeeze(x9m(:,1,xmc)),'-','linewidth',1.5)

ylim([-eulLim eulLim])

legend 'phi true' 'phi Ayaw Acc' 'phi Byaw Acc' 'theta Ayaw+yawnoise'
subplot(2,1,2)
plot(clock,thetanoise,'linewidth',1.5)
hold on
plot(clock,squeeze(x7m(:,2,xmc)),'-','linewidth',1.5)
plot(clock,squeeze(x8m(:,2,xmc)),'-.','linewidth',1.5,'color','g')
plot(clock,squeeze(x9m(:,2,xmc)),'-','linewidth',1.5)
legend 'theta true' 'theta Ayaw Acc' 'theta Byaw Acc' 'theta Ayaw+yawnoise'
ylim([-eulLim eulLim])


m_x7m = mean(x7m,3);
m_x8m = mean(x8m,3);
std_x7m =  sqrt(mean((x7m-m_x7m).^2,3));
std_x8m =  sqrt(mean((x8m-m_x8m).^2,3));

m_x9m = mean(x9m,3);
m_x10m = mean(x10m,3);
std_x9m =  sqrt(mean((x9m-m_x9m).^2,3));
std_x10m =  sqrt(mean((x10m-m_x10m).^2,3));


ind = 10;ind0 = 1;numsigma = 3;
clock2 = [clock(ind0:ind:end) fliplr(clock(ind0:ind:end))];
bound7 = [m_x7m(ind0:ind:end,:) - numsigma*std_x7m(ind0:ind:end,:);...
    m_x7m(end:-ind:ind0,:) + numsigma*std_x7m(end:-ind:ind0,:)];
bound8 = [m_x8m(ind0:ind:end,:) - numsigma*std_x8m(ind0:ind:end,:);...
    m_x8m(end:-ind:ind0,:) + numsigma*std_x8m(end:-ind:ind0,:)];

bound9 = [m_x9m(ind0:ind:end,:) - numsigma*std_x9m(ind0:ind:end,:);...
    m_x9m(end:-ind:ind0,:) + numsigma*std_x9m(end:-ind:ind0,:)];
bound10 = [m_x10m(ind0:ind:end,:) - numsigma*std_x10m(ind0:ind:end,:);...
    m_x10m(end:-ind:ind0,:) + numsigma*std_x10m(end:-ind:ind0,:)];

bound7 = [squeeze(x7m(ind0:ind:end,:,xmc)) - numsigma*std_x7m(ind0:ind:end,:);...
    squeeze(x7m(end:-ind:ind0,:,xmc)) + numsigma*std_x7m(end:-ind:ind0,:)];
bound8 = [squeeze(x8m(ind0:ind:end,:,xmc)) - numsigma*std_x8m(ind0:ind:end,:);...
    squeeze(x8m(end:-ind:ind0,:,xmc)) + numsigma*std_x8m(end:-ind:ind0,:)];

bound9 = [squeeze(x9m(ind0:ind:end,:,xmc)) - numsigma*std_x9m(ind0:ind:end,:);...
    squeeze(x9m(end:-ind:ind0,:,xmc)) + numsigma*std_x9m(end:-ind:ind0,:)];
bound10 = [squeeze(x10m(ind0:ind:end,:,xmc)) - numsigma*std_x10m(ind0:ind:end,:);...
    squeeze(x10m(end:-ind:ind0,:,xmc)) + numsigma*std_x10m(end:-ind:ind0,:)];

ind = 10;font = 12;
ff(2) = figure(2);clf
ff(2).Position = [horix+movex very-0.5*movey 1.5*figx 1.5*figy];
sgtitle 'Gyro bias estimates with noise-free yaw measurements'
h1 = subplot(2,1,1);
h1.FontSize = font;hold on
plot(clock(1:ind:end),bias_pp(1:ind:end),'linewidth',1.5)
hold on
% plot(clock, squeeze(bias_pm(:,1)))
plot(clock(1:ind:end),squeeze(x7m(1:ind:end,9,xmc)),'-','linewidth',2,'color','r')
plot(clock(1:ind:end),squeeze(x8m(1:ind:end,9,xmc)),'-.','linewidth',2,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound7(:,9)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound8(:,9)],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{\omega x}$', '$\hat{b}_{\omega x}$ augmented model', '$\hat{b}_{\omega x}$ reduced augmented model'},'Interpreter','latex','Location','ne','FontSize',font)
grid on
xlabel('Time [s]','FontSize',font)
ylabel('x gyro bias','FontSize',font)
ylim(0.05*[-beulLim beulLim]+mean(bias_pp)+0.005)
box on

h2 = subplot(2,1,2);
h2.FontSize = font;hold on
plot(clock(1:ind:end),bias_qq(1:ind:end),'linewidth',1.5)
hold on
% plot(clock, squeeze(bias_qm(:,1)))
plot(clock(1:ind:end),squeeze(x7m(1:ind:end,10,xmc)),'-','linewidth',2,'color','r')
plot(clock(1:ind:end),squeeze(x8m(1:ind:end,10,xmc)),'-.','linewidth',2,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound7(:,10)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound8(:,10)],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{\omega y}$', '$\hat{b}_{\omega y}$ augmented model', '$\hat{b}_{\omega y}$ reduced augmented model'},'Interpreter','latex','Location','ne','FontSize',font)
grid on
xlabel('Time [s]','FontSize',font)
ylabel('y gyro bias','FontSize',font)
ylim(0.05*[-beulLim beulLim]+mean(bias_qq)+0.005)
box on

%
ff(3) = figure(3);clf
ff(3).Position = [horix+2*movex very-0.5*movey 1.5*figx 1.5*figy];
sgtitle 'Gyro bias estimates with noisy yaw measurements'
h1 = subplot(2,1,1);h1.FontSize = font;hold on
plot(clock(1:ind:end),bias_pp(1:ind:end),'linewidth',1.5)
hold on
% plot(clock, squeeze(bias_pm(:,1)))
plot(clock(1:ind:end),squeeze(x9m(1:ind:end,9,xmc)),'-','linewidth',2,'color','r')
plot(clock(1:ind:end),squeeze(x10m(1:ind:end,9,xmc)),'-.','linewidth',2,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound9(:,9)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound10(:,9)],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{\omega x}$', '$\hat{b}_{\omega x}$ augmented model', '$\hat{b}_{\omega x}$ reduced augmented model'},'Interpreter','latex','Location','ne','FontSize',font)
grid on
xlabel('Time [s]','FontSize',font)
ylabel('x gyro bias','FontSize',font)
ylim(0.05*[-beulLim beulLim]++mean(bias_pp)+0.005)
box on

h2 = subplot(2,1,2);h2.FontSize = font;hold on
plot(clock(1:ind:end),bias_qq(1:ind:end),'linewidth',1.5)
hold on
% plot(clock, squeeze(bias_qm(:,1)))
plot(clock(1:ind:end),squeeze(x9m(1:ind:end,10,xmc)),'-','linewidth',2,'color','r')
plot(clock(1:ind:end),squeeze(x10m(1:ind:end,10,xmc)),'-.','linewidth',2,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound9(:,10)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound10(:,10)],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{\omega y}$', '$\hat{b}_{\omega y}$ augmented model', '$\hat{b}_{\omega y}$ reduced augmented model'},'Interpreter','latex','Location','ne','FontSize',font)
grid on
xlabel('Time [s]','FontSize',font)
ylabel('y gyro bias','FontSize',font)
ylim(0.05*[-beulLim beulLim]+mean(bias_qq)+0.005)
box on



ff(4) = figure(4);clf
ff(4).Position = [horix+2.5*movex very-0.5*movey 1.5*figx 1.5*figy];
sgtitle 'Accelerometer bias estimates with noise free yaw measurements'
h1 = subplot(2,1,1);
h1.FontSize = font;hold on
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock(1:ind:end),squeeze(x7m(1:ind:end,7,xmc)),'-','linewidth',2,'color','r')
plot(clock(1:ind:end),squeeze(x8m(1:ind:end,7,xmc)),'-.','linewidth',2,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound7(:,7)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound8(:,7)],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{ax}$', '$\hat{b}_{ax}$ augmented model', '$\hat{b}_{ax}$ reduced augmented model'},'Interpreter','latex','Location','ne')
grid on
xlabel 'Time [s]'
ylabel 'x acc bias'
ylim([-baccLim baccLim]+mean(bias_ax)+0.2)
box on

h2 = subplot(2,1,2);
h2.FontSize = font;hold on
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock(1:ind:end),squeeze(x7m(1:ind:end,8,xmc)),'-','linewidth',2,'color','r')
plot(clock(1:ind:end),squeeze(x8m(1:ind:end,8,xmc)),'-.','linewidth',2,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound7(:,8)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound8(:,8)],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{ay}$', '$\hat{b}_{ay}$ augmented model', '$\hat{b}_{ay}$ reduced augmented model'},'Interpreter','latex','Location','ne')
grid on
xlabel 'Time [s]'
ylabel 'y acc bias'
ylim([-baccLim baccLim]+mean(bias_ay)+0.2)
box on


ff(5) = figure(5);clf
ff(5).Position = [horix very-movey 1.5*figx 1.5*figy];
sgtitle 'Accelerometer bias estimates with noisy yaw measurements'
h1 = subplot(2,1,1);
h1.FontSize = font;hold on
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock(1:ind:end),squeeze(x9m(1:ind:end,7,xmc)),'-','linewidth',2,'color','r')
plot(clock(1:ind:end),squeeze(x10m(1:ind:end,7,xmc)),'-.','linewidth',2,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound9(:,7)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound10(:,7)],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{ax}$', '$\hat{b}_{ax}$ augmented model', '$\hat{b}_{ax}$ reduced augmented model'},'Interpreter','latex','Location','ne')
grid on
xlabel 'Time [s]'
ylabel 'x acc bias'
ylim([-baccLim baccLim]+mean(bias_ax)+0.2)
box on

h2 = subplot(2,1,2);
h2.FontSize = font;hold on
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock(1:ind:end),squeeze(x9m(1:ind:end,8,xmc)),'-','linewidth',2,'color','r')
plot(clock(1:ind:end),squeeze(x10m(1:ind:end,8,xmc)),'-.','linewidth',2,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound9(:,8)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound10(:,8)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$b_{ay}$', '$\hat{b}_{ay}$ augmented model', '$\hat{b}_{ay}$ reduced augmented model'},'Interpreter','latex','Location','ne')
grid on
xlabel 'Time [s]'
ylabel 'y acc bias'
ylim([-baccLim baccLim]+mean(bias_ay)+0.2)
box on

%
font = 14;
ff(6) = figure(6);clf
ff(6).Position = [horix+movex very-movey 2.5*figx 2.5*figy];
sgtitle 'State estimates with noise-free yaw measurements'
h1 = subplot(3,2,1);h1.FontSize = font;
hold on
plot(clock,stateData(:,4)-m_x7m(:,1),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,4)-m_x8m(:,1),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,4);stateData(end:-ind:ind0,4)]-bound7(:,1)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,4);stateData(end:-ind:ind0,4)]-bound8(:,1)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$\phi-\hat{\phi}$ augmented model','$\phi-\hat{\phi}$ reduced augmented model'},'interpreter','latex')
xlabel('Time [s]')
ylabel('$\phi-\hat{\phi}$','interpreter','latex')
ylim(0.25*[-eulLim eulLim])
xlim([0 120])
grid on
box on

h2 = subplot(3,2,2);h2.FontSize = font;
hold on
plot(clock,stateData(:,5)-m_x7m(:,2),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,5)-m_x8m(:,2),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,5);stateData(end:-ind:ind0,5)]-bound7(:,2)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,5);stateData(end:-ind:ind0,5)]-bound8(:,2)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$\theta-\hat{\theta}$ augmented model','$\theta-\hat{\theta}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$\theta-\hat{\theta}$','interpreter','latex')
ylim(0.25*[-eulLim eulLim])
xlim([0 120])
box on

h3 = subplot(3,2,3);h3.FontSize = font;hold on
plot(clock,stateData(:,13)-m_x7m(:,3),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,13)-m_x8m(:,3),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,13);stateData(end:-ind:ind0,13)]-bound7(:,3)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,13);stateData(end:-ind:ind0,13)]-bound8(:,3)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$u-\hat{u}$ augmented model','$u-\hat{u}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$u-\hat{u}$','interpreter','latex')
ylim([-eulLim eulLim])
xlim([0 120])
box on

h4 = subplot(3,2,4);h4.FontSize = font;hold on
plot(clock,stateData(:,14)-m_x7m(:,4),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,14)-m_x8m(:,4),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,14);stateData(end:-ind:ind0,14)]-bound7(:,4)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,14);stateData(end:-ind:ind0,14)]-bound8(:,4)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$v-\hat{v}$ augmented model','$v-\hat{v}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$v-\hat{v}$','interpreter','latex')
ylim([-eulLim eulLim])
xlim([0 120])
box on

h5 = subplot(3,2,5);h5.FontSize = font;hold on
plot(clock,stateData(:,1)-m_x7m(:,5),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,1)-m_x8m(:,5),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,1);stateData(end:-ind:ind0,1)]-bound7(:,5)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,1);stateData(end:-ind:ind0,1)]-bound8(:,5)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$x-\hat{x}$ augmented model','$x-\hat{x}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$x-\hat{x}$','interpreter','latex')
ylim([-eulLim eulLim])
xlim([0 120])
box on

h6 = subplot(3,2,6);h6.FontSize = font;hold on
plot(clock,stateData(:,2)-m_x7m(:,6),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,2)-m_x8m(:,6),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,2);stateData(end:-ind:ind0,2)]-bound7(:,6)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,2);stateData(end:-ind:ind0,2)]-bound8(:,6)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$y-\hat{y}$ augmented model','$y-\hat{y}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$y-\hat{y}$','interpreter','latex')
ylim([-eulLim eulLim])
xlim([0 120])
box on

%
ff(7) = figure(7);clf
ff(7).Position = [horix+2*movex very-movey 2.5*figx 2.5*figy];
sgtitle 'State estimates with noisy yaw measurements'
h1 = subplot(3,2,1);h1.FontSize = font;
hold on
plot(clock,stateData(:,4)-m_x9m(:,1),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,4)-m_x10m(:,1),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,4);stateData(end:-ind:ind0,4)]-bound9(:,1)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,4);stateData(end:-ind:ind0,4)]-bound10(:,1)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$\phi-\hat{\phi}$ augmented model','$\phi-\hat{\phi}$ reduced augmented model'},'interpreter','latex')
xlabel('Time [s]')
ylabel('$\phi-\hat{\phi}$','interpreter','latex')
ylim(0.25*[-eulLim eulLim])
xlim([0 120])
grid on
box on

h2 = subplot(3,2,2);h2.FontSize = font;
hold on
plot(clock,stateData(:,5)-m_x9m(:,2),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,5)-m_x10m(:,2),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,5);stateData(end:-ind:ind0,5)]-bound9(:,2)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,5);stateData(end:-ind:ind0,5)]-bound10(:,2)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$\theta-\hat{\theta}$ augmented model','$\theta-\hat{\theta}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$\theta-\hat{\theta}$','interpreter','latex')
ylim(0.25*[-eulLim eulLim])
xlim([0 120])
box on

h3 = subplot(3,2,3);h3.FontSize = font;hold on
plot(clock,stateData(:,13)-m_x9m(:,3),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,13)-m_x10m(:,3),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,13);stateData(end:-ind:ind0,13)]-bound9(:,3)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,13);stateData(end:-ind:ind0,13)]-bound10(:,3)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$u-\hat{u}$ augmented model','$u-\hat{u}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$u-\hat{u}$','interpreter','latex')
ylim([-eulLim eulLim])
xlim([0 120])
box on

h4 = subplot(3,2,4);h4.FontSize = font;hold on
plot(clock,stateData(:,14)-m_x9m(:,4),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,14)-m_x10m(:,4),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,14);stateData(end:-ind:ind0,14)]-bound9(:,4)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,14);stateData(end:-ind:ind0,14)]-bound10(:,4)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$v-\hat{v}$ augmented model','$v-\hat{v}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$v-\hat{v}$','interpreter','latex')
ylim([-eulLim eulLim])
xlim([0 120])
box on

h5 = subplot(3,2,5);h5.FontSize = font;hold on
plot(clock,stateData(:,1)-m_x9m(:,5),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,1)-m_x10m(:,5),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,1);stateData(end:-ind:ind0,1)]-bound9(:,5)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,1);stateData(end:-ind:ind0,1)]-bound10(:,5)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$x-\hat{x}$ augmented model','$x-\hat{x}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$x-\hat{x}$','interpreter','latex')
ylim([-eulLim eulLim])
xlim([0 120])
box on

h6 = subplot(3,2,6);h6.FontSize = font;hold on
plot(clock,stateData(:,2)-m_x9m(:,6),'-','linewidth',1.5,'color','r')
plot(clock,stateData(:,2)-m_x10m(:,6),'-.','linewidth',1.5,'color','g')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,2);stateData(end:-ind:ind0,2)]-bound9(:,6)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.3);
patch('Faces',(1:length(clock2)),'Vertices',[clock2' [stateData(ind0:ind:end,2);stateData(end:-ind:ind0,2)]-bound10(:,6)],'EdgeColor','none','FaceColor','green','LineWidth',1,'FaceAlpha',0.2);
legend({'$y-\hat{y}$ augmented model','$y-\hat{y}$ reduced augmented model'},'interpreter','latex')
grid on
xlabel('Time [s]')
ylabel('$y-\hat{y}$','interpreter','latex')
ylim([-eulLim eulLim])
xlim([0 120])
box on

%
tic
ff(8) = figure(8);clf
title('Sensor bias estimation','FontSize',12)
ff(8).Position = [horix+2*movex very-0.5*movey 2*figx 2*figy];
ff(1,1) = subplot(2,2,1);ff(1,1).FontSize = 12;hold on
plot(clock(1:ind:end),bias_pp(1:ind:end),'LineWidth',1.5)
hold on
plot(clock(1:ind:end),m_x10m(1:ind:end,9),'-.','linewidth',2,'color','r')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound8(:,9)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{\omega x}$', '$\hat{b}_{\omega x}$', '$\sigma \hat{b}_{\omega x}$'},'Interpreter','latex','Location','ne','FontSize',font)
grid on
xlabel('Time [s]','FontSize',font)
ylabel('x gyro bias','FontSize',font)
ylim(0.05*[-beulLim beulLim]+mean(bias_pp)+0.005)
xlim([0 120])
box on

ff(1,2) = subplot(2,2,2);ff(1,2).FontSize = font;hold on
plot(clock(1:ind:end),bias_qq(1:ind:end),'linewidth',1.5)
hold on
plot(clock(1:ind:end),m_x10m(1:ind:end,10),'-.','linewidth',2,'color','r')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound8(:,10)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{\omega y}$', '$\hat{b}_{\omega y}$', '$\sigma \hat{b}_{\omega y}$'},'Interpreter','latex','Location','ne','FontSize',font)
grid on
xlabel('Time [s]','FontSize',font)
ylabel('y gyro bias','FontSize',font)
ylim(0.05*[-beulLim beulLim]+mean(bias_qq)+0.005)
box on
xlim([0 120])

ff(2,1) = subplot(2,2,3);ff(2,1).FontSize = font;hold on
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,m_x10m(:,7),'-.','linewidth',2,'color','r')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound10(:,7)],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.2);
legend({'$b_{ax}$', '$\hat{b}_{ax}$', '$\sigma \hat{b}_{ax}$'},'Interpreter','latex','Location','ne')
grid on
xlabel('Time [s]','FontSize',font)
ylabel('x acc bias','FontSize',font)
ylim([-baccLim baccLim]+mean(bias_ax)+0.2)
box on
xlim([0 120])

ff(2,2) = subplot(2,2,4);ff(2,2).FontSize = font;hold on
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,m_x10m(:,8),'-.','linewidth',2,'color','r')
patch('Faces',(1:length(clock2)),'Vertices',[clock2' bound10(:,8)],'EdgeColor','none','FaceColor','red','LineWidth',1,'FaceAlpha',0.2);
legend({'$b_{ay}$', '$\hat{b}_{ay}$', '$\sigma \hat{b}_{ay}$'},'Interpreter','latex','Location','ne')
grid on
xlabel('Time [s]','FontSize',font)
ylabel('y acc bias','FontSize',font)
ylim([-baccLim baccLim]+mean(bias_ay)+0.2)
box on
xlim([0 120])






%%
% Now we consider the sensor-to-sensor model for mass change detection
% param.m = 0.5;
Gs = -tf([0 param.g(3)*param.lambda1/param.m],[1 param.lambda1/param.m 0]);
Gd = c2d(Gs,dt,'zoh');
[bhat,ahat] = tfdata(Gd,'v');
th = [ahat(2:end) 0 (bhat(2:end))]';

Nw = 2000;
lambdaW = ones(Nw,1);
for ii = 1:Nw
    lambdaW(ii,:) = 0.999^(Nw-ii);
end
lambdaW = sqrt(lambdaW);


q0 = IMU0w0psi(:,5);
ax0 = IMU0w0psi(:,1);
p0 = -IMU0w0psi(:,4);
ay0 = IMU0w0psi(:,2);

% Linear interpolation of the bias
bias_gyroout = x10m(:,9:10,xmc);
bias_accout = x10m(:,7:8,xmc);
xgyro = bias_gyroout;
xacc = bias_accout;

for ii = 1:N
    if mod(ii-1,index) == 0
        if ii<N-index+1
            rate = (bias_gyroout(ii+index,:)-bias_gyroout(ii,:))/index;
            xgyro(ii+1:ii+index-1,:) = xgyro(ii,:)+rate.*(1:index-1)';
            rate = (bias_accout(ii+index,:)-bias_accout(ii,:))/index;
            xacc(ii+1:ii+index-1,:) = xacc(ii,:)+rate.*(1:index-1)';
        end
    end
end


% With bias compensation
phatpsi = -IMU1w1psi(:,4) + (bias_p - xgyro(:,1)) ;
qhatpsi = IMU1w1psi(:,5) + (bias_q - xgyro(:,2));%xgyro(:,2)
axpsi = IMU1w1psi(:,1) + (bias_ax - xacc(:,1));
aypsi = IMU1w1psi(:,2) + (bias_ay - xacc(:,2));

% No compensation
phat = -IMU1w1psi(:,4) + 1*(bias_p - 0*xgyro(:,1)) ;
qhat = IMU1w1psi(:,5) + 1*(bias_q - 0*xgyro(:,2));
ax = IMU1w1psi(:,1) + 1*(bias_ax - 0*xacc(:,1));
ay = IMU1w1psi(:,2) + 1*(bias_ay - 0*xacc(:,2));

Zx0 = [-ax0(2:end-1) -ax0(1:end-2) q0(3:end) q0(2:end-1) q0(1:end-2)];
Phix = [-ax(2:end-1) -ax(1:end-2) qhat(3:end) qhat(2:end-1) qhat(1:end-2)];
Yx = ax(3:end);

Zy0 = [-ay0(2:end-1) -ay0(1:end-2) p0(3:end) p0(2:end-1) p0(1:end-2)];
Phiy = [-ay(2:end-1) -ay(1:end-2) phat(3:end) phat(2:end-1) phat(1:end-2)];
Yy = ay(3:end);

Zx0psi = [-ax0(2:end-1) -ax0(1:end-2) q0(3:end) q0(2:end-1) q0(1:end-2)];
Phixpsi = [-axpsi(2:end-1) -axpsi(1:end-2) qhatpsi(3:end) qhatpsi(2:end-1) qhatpsi(1:end-2)];
Yxpsi = axpsi(3:end);

Zy0psi = [-ay0(2:end-1) -ay0(1:end-2) p0(3:end) p0(2:end-1) p0(1:end-2)];
Phiypsi = [-aypsi(2:end-1) -aypsi(1:end-2) phatpsi(3:end) phatpsi(2:end-1) phatpsi(1:end-2)];
Yypsi = aypsi(3:end);

Nn = size(Zy0,1);
for ii = Nw:Nn
    Z0yii = [Zy0(ii-Nw+1:ii,:)];% [zeros(1,5);Zy0(ii-Nw+1:ii-1,:)] [zeros(2,5);Zy0(ii-Nw+1:ii-2,:)],...
    %          [zeros(3,5);Zy0(ii-Nw+1:ii-3,:)] [zeros(4,5);Zy0(ii-Nw+1:ii-4,:)] [zeros(5,5);Zy0(ii-Nw+1:ii-5,:)]];
    Phiyii = Phiy(ii-Nw+1:ii,:);
    Yyii = Yy(ii-Nw+1:ii,:);
    Eyii = Yyii-Phiyii*th;
    YY = norm(Z0yii'*(Eyii.*lambdaW));
    Jyii(ii,:) = YY;
    
    Z0xii = Zx0(ii-Nw+1:ii,:);
    Phixii = Phix(ii-Nw+1:ii,:);
    Yxii = Yx(ii-Nw+1:ii,:);
    Exii = Yxii-Phixii*th;
    XX = norm(Z0xii'*(Exii.*lambdaW));
    Jxii(ii,:) = XX;
    
    Z0yiipsi = [Zy0psi(ii-Nw+1:ii,:)];
    Phiyiipsi = Phiypsi(ii-Nw+1:ii,:);
    Yyiipsi = Yypsi(ii-Nw+1:ii,:);
    Eyiipsi = Yyiipsi-Phiyiipsi*th;
    YYpsi = norm(Z0yiipsi'*(Eyiipsi.*lambdaW));
    Jyiipsi(ii,:) = YYpsi;
    
    Z0xiipsi = Zx0psi(ii-Nw+1:ii,:);
    Phixiipsi = Phixpsi(ii-Nw+1:ii,:);
    Yxiipsi = Yxpsi(ii-Nw+1:ii,:);
    Exiipsi = Yxiipsi-Phixiipsi*th;
    XXpsi = norm(Z0xiipsi'*(Exiipsi.*lambdaW));
    Jxiipsi(ii,:) = XXpsi;
end
%
ind = 10;
horix = -1900 ;very = 700;movex = 450; movey = 450; figx = 420; figy = 320;%560 420
h9 = figure;clf
h9.Position = [-1257         384        1.75*figx         figy];
h99 = axes;h99.FontSize = font;hold on
ff(1) = subplot(1,2,1);ff(1).FontSize = font;hold on
plot(clock(1,1:ind:Nn),Jxii(1:ind:end),'linewidth',2)
hold on
plot(clock(1,1:ind:Nn),Jyii(1:ind:end),'--r','linewidth',3)
xlabel 'Time [s]'
ylabel 'IV cost function J_{IV}'
legend 'Roll model' 'Pitch Model'
grid on
title('Without bias compensation','FontWeight','Normal')
ylim([0 0.1])
xlim([0 120])
box on

ff(2) = subplot(1,2,2);ff(2).FontSize = font;hold on
plot(clock(1,1:ind:Nn),Jxiipsi(1:ind:end),'linewidth',2)
hold on
plot(clock(1,1:ind:Nn),Jyiipsi(1:ind:end),'--r','linewidth',3)
xlabel 'Time [s]'
ylabel 'IV cost function J_{IV}'
legend 'Roll model' 'Pitch Model'
grid on
title('With bias compensation','FontWeight','Normal')
ylim([0 0.1])
xlim([0 120])
box on

%%
% Now we consider full model using a as output and input
% x = [phi theta u v x y biasphi biastheta biasx biasy windx windy]
% Tuning coefficient
Qeul = 1e-10*[1 1];
Qv = 1e-10*[1 1];
Qx = 1e-10*[1 1];
Qbiasg = 5*1e-5*[1 1];
Qbiasa = 5*1e-6*[1 1];
Qwind = 1e-4*[1 1];
index = 200;
clear x9m y9m x10m y10m P9m P10m P9eig P10eig
F9 = @(psit, r) eye(12)+param.dt*[  0 r, 0 0, 0 0, -1 0, 0 0, 0 0;...
    -r 0, 0 0, 0 0, 0 -1, 0 0, 0 0;...
    0 param.g(3), -param.lambda1/param.m r, 0 0, 0 0, 0 0, -param.lambda1/param.m*[cos(psit) sin(psit)] ;...
    -param.g(3) 0, -r -param.lambda1/param.m, 0 0, 0 0, 0 0, -param.lambda1/param.m*[-sin(psit) cos(psit)];...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0];

% F10 is nonlinear taking the cross production
% x = [phi theta u v x y biasphi biastheta biasx biasy]
F10 = @(psit, r) eye(10)+param.dt*[  0 r, 0 0, 0 0, -1 0, 0 0;...
    -r 0, 0 0, 0 0, 0 -1, 0 0;...
    0 param.g(3), 0 r, 0 0, 0 0, 1 0;...
    -param.g(3) 0, -r 0, 0 0, 0 0, 0 1;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];
H91 = @(psit) [0 0, 0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0, 0 0];
H92 = @(psit) [0 0, param.lambda1/param.m 0, 0 0, 0 0, 1 0, param.lambda1/param.m*[cos(psit) sin(psit)];...
    0 0, 0 param.lambda1/param.m, 0 0, 0 0, 0 1, param.lambda1/param.m*[-sin(psit) cos(psit)]];
H101 = @(psit) [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0];
H102 = [];

B9 = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0, 0 0]';
B10 = param.dt*[1 0, 0 0, 0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, -1 0, 0 0, 0 0, 0 0;...
    0 0, 0 -1, 0 0, 0 0, 0 0]';
G9 = [1 0, 0 0,0 0, 0 0,0 0,0 0;...
    0 1, 0 0, 0 0, 0 0,0 0,0 0;...
    0 0, param.dt 0, 0 0, 0 0,0 0,0 0;...
    0 0, 0 param.dt, 0 0, 0 0,0 0,0 0;...
    0 0, param.dt^2/2 0, 0 0, 0 0,0 0,0 0;...
    0 0, 0 param.dt^2/2, 0 0, 0 0,0 0,0 0;...
    0 0, 0 0,0 0, 1 0, 0 0,0 0;...
    0 0, 0 0,0 0, 0 1, 0 0,0 0;...
    0 0, 0 0,0 0, 0 0, 1 0,0 0;...
    0 0, 0 0,0 0, 0 0, 0 1,0 0;...
    0 0, 0 0,0 0, 0 0, 0 0,1 0;...
    0 0, 0 0,0 0, 0 0, 0 0,0 1];
G10 = [1 0, 0 0,0 0, 0 0, 0 0;...
    0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, param.dt 0, 0 0, 0 0, 0 0 ;...
    0 0, 0 param.dt, 0 0, 0 0, 0 0;...
    0 0, 0 0, param.dt^2/2 0, 0 0, 0 0;...
    0 0, 0 0, 0 param.dt^2/2, 0 0, 0 0;...
    0 0, 0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 0, 0 1, 0 0;...
    0 0, 0 0, 0 0, 0 0, 1 0;...
    0 0, 0 0, 0 0, 0 0, 0 1];



R91 = diag([Rpos Rpos Rspeed Rspeed]);
R92 = diag([Racc Racc]);
R9 = blkdiag(R91,R92);
R101 = R91;
R102 = diag([]);
R10 = blkdiag(R101,R102);
Q9 = diag([Qeul, Qv, Qx, Qbiasg, Qbiasa, Qwind]);
Q10 = diag([100*Qeul, Qv, Qx, 10*Qbiasg, 10*Qbiasa]);
% Stationary Kalman
% [K7bar,Pp7,Pf7] = dlqe(F7,G7,H7,Q7,R7);

x9 = zeros(12,1); x10 = zeros(10,1);
P9 = 1000*eye(12); P10 = 1000*eye(10);
for ii = 1:length(clock)
    % Data
    u9 = [pnoise(ii)+bias_p(ii) qnoise(ii)+bias_q(ii)];
    y9 = [xnoise(ii) ynoise(ii) xdotnoise(ii) ydotnoise(ii),...
        axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii)];
    u10 = [pnoise(ii)+bias_p(ii) qnoise(ii)+bias_q(ii) axnoise(ii)+bias_ax(ii) aynoise(ii)+bias_ay(ii)];
    y10 = [xnoise(ii) ynoise(ii),...
        xdotnoise(ii) ydotnoise(ii)];
    
    % Prediction step
    F9ii = F9(psinoise(ii),psidot(ii));
    x9 = F9ii*x9 + B9*u9';
    P9 = F9ii*P9*F9ii' + G9*Q9*G9';
    % Prediction step
    F10ii = F10(psinoise(ii),psidot(ii));
    x10 = F10ii*x10 + B10*u10';
    P10 = F10ii*P10*F10ii' + G10*Q10*G10';
    
    H91ii = H91(psinoise(ii));
    H92ii = H92(psinoise(ii));
    H9 = [H91ii;H92ii];
    
    H101ii = H101(psinoise(ii));
    H102ii = H102;
    H10 = [H101ii;H102ii];
    rr9(ii,:) = rank(obsv(F9ii,H9));
    cond9(ii,:) = cond(obsv(F9ii,H9));
    rr10(ii,:) = rank(obsv(F10ii,H10));
    cond10(ii,:) = cond(obsv(F10ii,H10));
    if mod(ii-1,index) == 0
        % Correction step
        S9 = H9*P9*H9'+R9;
        K9 = P9*H9'*inv(S9);
        x9 = x9 + K9*(y9(1:end)' - H9*x9);
        P9 = P9 - K9*H9*P9;
        
        S10 = H10*P10*H10'+R10;
        K10 = P10*H10'*inv(S10);
        x10 = x10 + K10*(y10(1:end)' - H10*x10);
        P10 = P10 - K10*H10*P10;
    else
        % Correction step
        S9 = H92ii*P9*H92ii'+R92;
        K9 = P9*H92ii'*inv(S9);
        x9 = x9 + K9*(y9(5:end)' - H92ii*x9);
        P9 = P9 - K9*H92ii*P9;
        
        %        S10 = H102ii*P10*H102ii'+R102;
        %        K10 = P10*H102ii'*inv(S10);
        %        x10 = x10 + K10*(y10(5:end)' - H102ii*x10);
        %        P10 = P10 - K10*H102ii*P10;
    end
    
    P9m(ii,:) = diag(P9);
    P9eig(ii,:) = eig(P9);
    P10m(ii,:) = diag(P10);
    P10eig(ii,:) = eig(P10);
    x9m(ii,:) = x9;
    y9m(ii,:) = H9*x9;
    
    x10m(ii,:) = x10;
    y10m(ii,:) = H10*x10;
end
disp('Done compute gyro and acc biases')

%
close all
baccLim = 0.6;windxLim = 4;windyLim = 2;uLim = 3;vLim = 3;

ff(5) = figure(5);clf
ff(5).Position = [horix very-movey figx figy];
subplot(2,1,1)
plot(clock,windData(:,1),'linewidth',1.5)
hold on
plot(clock,x9m(:,11),'-','linewidth',1.5)
legend 'windx true' 'windx aout'
ylim([-windxLim windxLim])

subplot(2,1,2)
plot(clock,windData(:,2),'linewidth',1.5)
hold on
plot(clock,x9m(:,12),'-','linewidth',1.5)
legend 'windy true'  'windy aout'
ylim([-windyLim windyLim])

ff(6) = figure(6);clf
ff(6).Position = [horix+movex very-movey figx figy];
subplot(2,1,1)
plot(clock,abs(stateData(:,4)-x9m(:,1)),'linewidth',1.5)
hold on
plot(clock,abs(stateData(:,4)-x10m(:,1)),'linewidth',1.5,'color','g')
legend 'Diff phi acc ain' 'Diff phi acc ain+non'
ylim([-eulLim eulLim])

subplot(2,1,2)
plot(clock,abs(stateData(:,5)-x9m(:,2)),'linewidth',1.5)
hold on
plot(clock,abs(stateData(:,5)-x10m(:,2)),'linewidth',1.5,'color','g')
legend 'Diff theta acc ain' 'Diff theta acc ain+non'
ylim([-eulLim eulLim])


ff(7) = figure(7);clf
ff(7).Position = [horix+2*movex very-movey figx figy];
subplot(2,1,1)
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,x9m(:,9),'-.','linewidth',1.5)
plot(clock,x10m(:,9),'--','linewidth',1.5,'color','g')
legend 'bias ax true' 'bias ax ain' 'bias ax ain+non'
ylim([-baccLim baccLim])


subplot(2,1,2)
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,x9m(:,10),'-.','linewidth',1.5)
plot(clock,x10m(:,10),'--','linewidth',1.5,'color','g')
legend 'bias ay true' 'bias ay ain' 'bias ay ain+non'
ylim([-baccLim baccLim])

ff(8) = figure(8);clf
ff(8).Position = [horix+3*movex very-movey figx figy];
subplot(2,1,1)
plot(clock,bias_p,'linewidth',1.5)
hold on
plot(clock,x9m(:,7),'linewidth',1.5)
plot(clock,x10m(:,7),'linewidth',1.5,'color','g')
ylim([-eulLim eulLim])
subplot(2,1,2)
plot(clock,bias_q,'linewidth',1.5)
hold on
plot(clock,x9m(:,8),'linewidth',1.5)
plot(clock,x10m(:,8),'linewidth',1.5,'color','g')
ylim([-eulLim eulLim])

ccii = 0.1*[ones(1,500) -ones(1,500) -ones(1,500) ones(1,500)];

%
ff(5) = figure(5);clf % bias ax
ff(5).Position = [horix very-movey figx figy];
subplot(2,1,1)
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,x3m(:,7),'--','linewidth',1.5)
plot(clock,x4m(:,7),'-.','linewidth',1.5)
ylim([-baccLim baccLim])
legend 'b-a_x' 'b-a_x Aout Acc' 'b-a_x Ain Acc'
xlabel 'Time [s]'
grid on


subplot(2,1,2)
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,x3m(:,8),'--','linewidth',1.5)
plot(clock,x4m(:,8),'-.','linewidth',1.5)
ylim([-baccLim baccLim])
legend 'b-a_y' 'b-a_y Aout Acc' 'b-a_y Aout Acc'
xlabel 'Time [s]'
grid on


ff(6) = figure(6);clf % bias ay
ff(6).Position = [horix+movex very-movey figx figy];
subplot(2,1,1)
plot(clock,bias_ax,'linewidth',1.5)
hold on
plot(clock,x5m(:,9),'-','linewidth',1.5)
plot(clock,x6Smooth(:,9),'-.','linewidth',1.5)
ylim([-baccLim baccLim])
legend 'b-a_x' 'b-a_x Aout Gyro+Acc' 'b-a_x Ain Gyro+Acc'
xlabel 'Time [s]'
grid on


subplot(2,1,2)
plot(clock,bias_ay,'linewidth',1.5)
hold on
plot(clock,x5m(:,10),'-','linewidth',1.5)
plot(clock,x6Smooth(:,10),'-.','linewidth',1.5)
ylim([-baccLim baccLim])
legend 'b-a_y' 'b-a_y Aout Gyro+Acc' 'b-a_y Aout Gyro+Acc'
xlabel 'Time [s]'
grid on


%
factor = 1;
ff(7) = figure(7);clf % bias p
ff(7).Position = [horix+2*movex very-movey figx figy];
h7 = axes;h7.FontSize = 14;hold on
subplot(2,1,1)
plot(clock(1:factor:end),bias_p(1:factor:end),'linewidth',1.5)
hold on
plot(clock(1:factor:end),x5m(1:factor:end,7),'--','linewidth',1.5)
plot(clock(1:factor:end),x6m(1:factor:end,7),'-.','linewidth',1.5)
ylim([-beulLim beulLim])
xlabel 'Time [s]'
ylabel '\beta_\phi'
legend 'b-phi true' 'b-phi Aout GYRO+Acc' 'b-phi Ain GYRO+Acc'
grid on

subplot(2,1,2)
plot(clock(1:factor:end),bias_q(1:factor:end),'linewidth',1.5)
hold on
plot(clock(1:factor:end),x5m(1:factor:end,8),'--','linewidth',1.5)
plot(clock(1:factor:end),x6m(1:factor:end,8),'-.','linewidth',1.5)
ylim([-beulLim beulLim])
xlabel 'Time [s]'
ylabel '\beta_\theta'
legend 'b-theta true' 'b-theta Aout GYRO+Acc' 'b-theta Ain GYRO+Acc'
grid on


ff(8) = figure(8);clf % bias q
ff(8).Position = [horix+3*movex very-movey figx figy];
h8 = axes;h8.FontSize = 14;hold on
hold on




%% Consider the no wind model
% Model 1: x = [phi,theta,u,v,x,y,betaphi,betatheta]; y = [x,y,xdot,ydot,ax,ay]
% model 2: x = [phi,theta,u,v,x,y,betaphi,betatheta]; y = [x,y,xdot,ydot]
% model 3: x = [phi,theta,u,v,x,y,betax,betay]; y = [x,y,xdot,ydot,ax,ay]
% model 4: x = [phi,theta,u,v,x,y,betax,betay]; y = [x,y,xdot,ydot]
% model 5: x = [phi,theta,u,v,x,y,betaphi,betatheta,betax,betay]; y = [x,y,xdot,ydot,ax,ay]
% model 6: x = [phi,theta,u,v,x,y,betaphi,betatheta,betax,betay]; y = [x,y,xdot,ydot]

% And similar for the yawrate model

% GYRO
F11 = eye(8) + param.dt*[0 0, 0 0, 0 0, -1 0;...
    0 0, 0 0, 0 0, 0 -1;...
    0 -param.g(3), -param.lambda1/param.m 0, 0 0, 0 0;...
    param.g(3) 0, 0 -param.lambda1/param.m, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H11 = [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0;...
    0 0, param.lambda1/param.m 0, 0 0, 0 0;...
    0 0, 0 param.lambda1/param.m, 0 0, 0 0];

F12 = eye(8) + param.dt*[0 0, 0 0, 0 0, -1 0;...
    0 0, 0 0, 0 0, 0 -1;...
    0 -param.g(3), 0 0, 0 0, 0 0;...
    param.g(3) 0, 0 0, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H12 = [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0];


% ACC
F13 = eye(8) + param.dt*[0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 -param.g(3), -param.lambda1/param.m 0, 0 0, 0 0;...
    param.g(3) 0, 0 -param.lambda1/param.m, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H13 = [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0;...
    0 0, param.lambda1/param.m 0, 0 0, 1 0;...
    0 0, 0 param.lambda1/param.m, 0 0, 0 1];

F14 = eye(8) + param.dt*[0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 -param.g(3), 0 0, 0 0, -1 0;...
    param.g(3) 0, 0 0, 0 0, 0 -1;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H14 = [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, 1 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0];

% GYRO+ACC
% ACC
F15 = eye(10) + param.dt*[0 0, 0 0, 0 0, -1 0, 0 0;...
    0 0, 0 0, 0 0, 0 -1, 0 0;...
    0 -param.g(3), -param.lambda1/param.m 0, 0 0, 0 0, 0 0;...
    param.g(3) 0, 0 -param.lambda1/param.m, 0 0, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];
H15 = [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, param.lambda1/param.m 0, 0 0, 0 0, 1 0;...
    0 0, 0 param.lambda1/param.m, 0 0, 0 0, 0 1];

F16 = eye(10) + param.dt*[0 0, 0 0, 0 0, -1 0, 0 0;...
    0 0, 0 0, 0 0, 0 -1, 0 0;...
    0 -param.g(3), 0 0, 0 0, 0 0, 1 0;...
    param.g(3) 0, 0 0, 0 0, 0 0, 0 1;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];
H16 = [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0];


ii = randi(N);

% Now consider yaw rate
psi = stateNon0w1psi(:,6);
psidot = stateNon0w1psi(:,12);

psit = psi(ii);
rt = psidot(ii);
% GYRO
F17 = eye(8) + param.dt*[0 rt, 0 0, 0 0, -1 0;...
    -rt 0, 0 0, 0 0, 0 -1;...
    0 -param.g(3), -param.lambda1/param.m rt, 0 0, 0 0;...
    param.g(3) 0, -rt -param.lambda1/param.m, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H17 = [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, param.lambda1/param.m 0, 0 0, 0 0;...
    0 0, 0 param.lambda1/param.m, 0 0, 0 0];

F18 = eye(8) + param.dt*[0 rt, 0 0, 0 0, -1 0;...
    -rt 0, 0 0, 0 0, 0 -1;...
    0 -param.g(3), 0 rt, 0 0, 0 0;...
    param.g(3) 0, -rt 0, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H18 = [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0];


% ACC
F19 = eye(8) + param.dt*[0 rt, 0 0, 0 0, 0 0;...
    -rt 0, 0 0, 0 0, 0 0;...
    0 -param.g(3), -param.lambda1/param.m rt, 0 0, 0 0;...
    param.g(3) 0, -rt -param.lambda1/param.m, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H19 = [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, param.lambda1/param.m 0, 0 0, 1 0;...
    0 0, 0 param.lambda1/param.m, 0 0, 0 1];

F20 = eye(8) + param.dt*[0 rt, 0 0, 0 0, 0 0;...
    -rt 0, 0 0, 0 0, 0 0;...
    0 -param.g(3), 0 rt, 0 0, 1 0;...
    param.g(3) 0, -rt 0, 0 0, 0 1;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0];
H20 = [0 0, 0 0, 1 0, 0 0;...
    0 0, 0 0, 0 1, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0];

% GYRO+ACC
% ACC
F21 = eye(10) + param.dt*[0 rt, 0 0, 0 0, -1 0, 0 0;...
    -rt 0, 0 0, 0 0, 0 -1, 0 0;...
    0 -param.g(3), -param.lambda1/param.m rt, 0 0, 0 0, 0 0;...
    param.g(3) 0, -rt -param.lambda1/param.m, 0 0, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];
H21 = [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, param.lambda1/param.m 0, 0 0, 0 0, 1 0;...
    0 0, 0 param.lambda1/param.m, 0 0, 0 0, 0 1];

F22 = eye(10) + param.dt*[0 rt, 0 0, 0 0, -1 0, 0 0;...
    -rt 0, 0 0, 0 0, 0 -1, 0 0;...
    0 -param.g(3), 0 rt, 0 0, 0 0, 1 0;...
    param.g(3) 0, -rt 0, 0 0, 0 0, 0 1;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0];
H22 = [0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0];


% Using constant acceleration model
% x23 = [phi theta au av u v x y beta_phi beta_theta beta_x beta_y]
F23 = eye(12) + param.dt*[0 0, 0 0, 0 0, 0 0, -1 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 -1, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 1 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 1, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0];
%y23 = [x y xdot ydot ax ay] ax = -au-g*theta ay = g*phi-av
H23 = [0 0, 0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, 1 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 1, 0 0, 0 0, 0 0;...
    0 -param.g(3), -1 0, 0 0, 0 0, 0 0, -1 0;...
    param.g(3) 0, 0 -1, 0 0, 0 0, 0 0, 0 -1];


% Using constant acceleration model with psi change
% x23 = [phi theta au av u v x y beta_phi beta_theta beta_x beta_y]
F24 = eye(12) + param.dt*[0 rt, 0 0, 0 0, 0 0, -1 0, 0 0;...
    -rt 0, 0 0, 0 0, 0 0, 0 -1, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 1 0, 0 rt, 0 0, 0 0, 0 0;...
    0 0, 0 1, -rt 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 0, 0 0, 0 0];
%y23 = [x y xdot ydot ax ay] ax = -au-g*theta ay = g*phi-av
H24 = [0 0, 0 0, 0 0, 1 0, 0 0, 0 0;...
    0 0, 0 0, 0 0, 0 1, 0 0, 0 0;...
    0 0, 0 0, cos(psit) -sin(psit), 0 0, 0 0, 0 0;...
    0 0, 0 0, sin(psit) cos(psit), 0 0, 0 0, 0 0;...
    0 -param.g(3), -1 0, 0 0, 0 0, 0 0, -1 0;...
    param.g(3) 0, 0 -1, 0 0, 0 0, 0 0, 0 -1];

[rank(obsv(F11,H11)),rank(obsv(F12,H12)),rank(obsv(F13,H13)),rank(obsv(F14,H14)),...
    rank(obsv(F15,H15)),rank(obsv(F16,H16));...
    rank(obsv(F17,H17)),rank(obsv(F18,H18)),rank(obsv(F19,H19)),rank(obsv(F20,H20)),...
    rank(obsv(F21,H21)),rank(obsv(F22,H22))]









%%
% CORRECT SIMULATION, NOW SHOW SOME FIGURES
figx = 420;figy = 320;
font  = 12;
f1 = figure(1);clf
f1.Position = [-1850         350        1.5*figx         1.5*figy];
h1 = axes
plot3(pref(:,1),pref(:,2),pref(:,3),'-b','linewidth',2);
hold on
plot3(stateNon0w0psi(:,1),stateNon0w0psi(:,2),stateNon0w0psi(:,3),'--r','linewidth',3);
% hold on
% for ii = 1:1500:N
%     p_0ii = [stateNon0w0psi(ii,1),stateNon0w0psi(ii,2),stateNon0w0psi(ii,3)];
%     psiii = stateNon0w1psi(ii,6)
%     p_1ii = p_0ii + 5*[cos(psiii),sin(psiii),0];
%     plot3([p_0ii(1) p_1ii(1)],[p_0ii(2) p_1ii(2)],[p_0ii(3) p_1ii(3)],'g','linewidth',2);
% end

h1.FontSize = font;
hold off
xlabel 'X[m]'
ylabel 'Y[m]'
zlabel 'Z[m]'
legend 'Desired positions' 'Quadcopter positions'
zlim([7 9])
h1.XTick = [0 5 10];
h1.YTick = [0 5 10];
h1.ZTick = [7 8 9];
grid on
box on

f2 = figure(2);clf
f2.Position = [-1150 50 1.5*figx 1.5*figy];
h21 = subplot(3,1,1);
h21.FontSize = font;
box on
hold on
plot(clock,pref(:,1),'-b','linewidth',2)
plot(clock,stateNon1w1psi(:,1),'-.r','linewidth',3)
% plot(clock,stateNonSim1w0psi(:,1),'-g','linewidth',2)
ylim([-20 20])
h21.YTick = [-15 0 15];
ylabel 'X[m]'
xlabel 'Times[s]'
legend('Desired x position','True x position','Orientation','horizontal','Location','southeast')

h22 = subplot(3,1,2);
h22.FontSize = font;
hold on
plot(clock,pref(:,2),'-b','linewidth',2)
plot(clock,stateNon1w1psi(:,2),'-.r','linewidth',3)
% plot(clock,stateNonSim1w0psi(:,2),'-g','linewidth',2)
ylim([-20 20])
h22.YTick = [-15 0 15];
ylabel 'Y[m]'
xlabel 'Times[s]'
legend('Desired y position','True y position','Orientation','horizontal','Location','southeast')
box on

h23 = subplot(3,1,3);
h23.FontSize = font;
hold on
plot(clock,pref(:,3),'-b','linewidth',2)

plot(clock,stateNon1w1psi(:,3),'-.r','linewidth',3)
% plot(clock,stateNonSim1w0psi(:,3),'-g','linewidth',2)
ylim([5 11])
h23.YTick = [6 8 10];
ylabel 'Z[m]'
xlabel 'Times[s]'
legend('Desired z position','True z position','Orientation','horizontal','Location','southeast')
box on

%%
f3 = figure(3);clf
f3.Position = [-1800 300 1100 750];
% f3.Title =  'Residual of roll rates from true model and simulated model'
h31 = subplot(2,2,1);
h31.FontSize = 14;
hold on
plot(clock,IMU1w0psi(:,4)-IMUSim1w0psi(:,4),'linewidth',1.5)
title('$p_m-\mathring{p}$','Interpreter','latex','FontSize',18)
xlim([0 120])
ylabel '[rad/s]'
xlabel 'Times[s]'



h32 = subplot(2,2,2);
h32.FontSize = 14;
hold on
plot(clock,IMU1w0psi(:,5)-IMUSim1w0psi(:,5),'linewidth',1.5)
ylabel '[rad/s]'
xlabel 'Times[s]'
title('$q_m-\mathring{q}$','Interpreter','latex','FontSize',18)
xlim([0 120])

h33 = subplot(2,2,3);
h33.FontSize = 14;
hold on
plot(clock,IMU1w0psi(:,1)-IMUSim1w0psi(:,1),'-b','linewidth',1.5)
ylabel '[m/s^2]'
xlabel 'Times[s]'
title('$a_x-\mathring{a}_x$','Interpreter','latex','FontSize',18)
xlim([0 120])

h34 = subplot(2,2,4);
h34.FontSize = 14;
hold on
plot(clock,IMU1w0psi(:,2)-IMUSim1w0psi(:,2),'-b','linewidth',1.5)
plot(clock,abs(IMU1w0psi(:,2)-IMUSim1w0psi(:,2)))
ylabel '[m/s^2]'
xlabel 'Times[s]'
title('$a_y-\mathring{a}_{y}$','Interpreter','latex','FontSize',18)
xlim([0 120])


%%
f1 = figure(1);clf
f1.Position = [-1750 150 1100 750];
subplot(3,1,1)
plot(TxL,VxL,'linewidth',2)
hold on
plot(clock,stateNon1w1psi(:,7),'--','linewidth',2)
legend('x- desired velocity','x- true velocity','fontsize',14)
title('The longitudinal velocity','fontsize',14)
xlabel('Time[s]','fontsize',14)
ylabel('V_x[m/s]','fontsize',14)

subplot(3,1,2)
plot(TxL,VyL,'linewidth',2)
hold on
plot(clock,stateNon1w1psi(:,8),'--','linewidth',2)
legend('y- desired velocity','y- true velocity','fontsize',14)
title('The lateral velocity','fontsize',14)
xlabel('Time[s]','fontsize',14)
ylabel('V_y[m/s]','fontsize',14)

subplot(3,1,3)
plot(clock,psiref1w1psi,'-','linewidth',2)
hold on
plot(clock,stateNon1w1psi(:,6),'--','linewidth',2)
legend('\psi-desired angle','\psi-true angle','fontsize',14)
title('Course angle','fontsize',14)
xlabel('Time[s]','fontsize',14)
ylabel('\psi[rad]','fontsize',14)
%% PLOT WIND VELOCITY
f2 = figure(2);clf
f2.Position = [250 300 1.5*figx 1.5*figy];
f3a = axes;f3a.FontSize = 18;hold on
h1 = subplot(3,1,1);
h1.FontSize = font;
hold on
plot(clock,wind1w1psi(:,1),'-b','linewidth',2)
xlabel 'Time[s]'
ylabel 'u_w [m/s]'
legend 'in X- Direction'
h1.YLim = [0.8 1.2];
% h3.YTick = [-0.2 0 0.2];
box on

h2 = subplot(3,1,2);
h2.FontSize = font;
hold on
plot(clock,wind1w1psi(:,2),'-b','linewidth',2)
xlabel 'Time[s]'
ylabel 'v_w [m/s]'
legend 'in Y- Direction'
h2.YLim = [-1.2 -0.8];
% h2.YTick = [-1.2 -0.8];
box on

h3 = subplot(3,1,3);
h3.FontSize = font;
hold on
plot(clock,wind1w1psi(:,3),'-b','linewidth',2)
xlabel 'Time[s]'
ylabel 'w_w [m/s]'
legend 'in Z- Direction'
h3.YLim = [-0.2 0.2];
% h3.YTick = [-0.2 0.2];
box on
%% PLOT EULER ANGLE
ind = 100;
f3 = figure(3);clf % Plot the course angle
f3.Position = [250 300 1.5*figx 1.5*figy];
f3a = axes;f3a.FontSize = font;hold on
plot(clock(1:ind:end),psiref1w1psi(1:ind:end)/pi*180,'-','linewidth',2)
hold on

plot(clock(1:ind:end),stateNon1w1psi(1:ind:end,6)/pi*180,'-.r','linewidth',2.5)
hold off
legend 'Command course angle' 'Course angle'
xlabel 'Time[s]'
ylabel 'Course angle in degree'
ylim([-50 50])
box on

%% PLOT THE ORIENTATION ESTIMATION
f4 = figure(4);clf
f4.Position = [250 300 1100 750];
subplot(2,2,1)
plot(clock, phiEst(:,1),'linewidth',2)
hold on
plot(clock, phiEst(:,2),'--','linewidth',2)
plot(clock, phiEst(:,3),'-.','linewidth',2)
plot(clock, phiEst(:,4),'--','linewidth',2)
% ylim([-0.2 0.2])
xlim([0 120])
xlabel('Time[s]','FontSize',18)
ylabel({'$\hat{\phi}$'},'Interpreter','latex','FontSize',18)
title('Estimation of the roll angle','FontSize',18)
legend 'True' 'Complementary' 'From high-pass filtered gyro' 'From low-pass filtered acc'
grid on

subplot(2,2,3)
plot(clock, phiEst(:,1)-phiEst(:,2),'linewidth',2)
ylim([-0.05 0.05])
xlim([0 120])
xlabel('Time[s]','FontSize',18)
ylabel({'$\Delta\hat{\phi}$'},'Interpreter','latex','FontSize',18)
title('Error of the roll angle estimation','FontSize',18)
grid on

subplot(2,2,2)
plot(clock, thetaEst(:,1),'linewidth',2)
hold on
plot(clock, thetaEst(:,2),'--','linewidth',2)
plot(clock, thetaEst(:,3),'-.','linewidth',2)
plot(clock, thetaEst(:,4),'--','linewidth',2)
% ylim([-0.2 0.2])
xlim([0 120])
xlabel('Time[s]','FontSize',18)
ylabel({'$\hat{\theta}$'},'Interpreter','latex','FontSize',18)
title('Estimation of the pitch angle','FontSize',18)
legend 'True' 'Complementary' 'From high-pass filtered gyro' 'From low-pass filtered acc'
grid on


subplot(2,2,4)
plot(clock, thetaEst(:,1)-thetaEst(:,2),'linewidth',2)
ylim([-0.05 0.05])
xlim([0 120])
xlabel('Time[s]','FontSize',18)
ylabel({'$\Delta\hat{\theta}$'},'Interpreter','latex','FontSize',18)
title('Error of the pitch angle estimation','FontSize',18)
grid on

%% PLOT XY estimated and phi theta estimated
f5 = figure(5);clf
f5.Position = figPos;
subplot(2,1,1)
plot(clock, stateNon1w0psi(:,1),'b','linewidth',2)
hold on
plot(clock, stateNon1w0psi(:,2),'r','linewidth',2)
plot(clock, stateLin1w0psi(:,1),'--b','linewidth',2)
plot(clock, stateLin1w0psi(:,2),'--r','linewidth',2)
hold off
xlabel 'Time[s]'
ylabel 'm'
legend 'X true' 'Y true' 'X simulate' 'Y simulate'

subplot(2,1,2)
plot(clock, stateNon1w0psi(:,4),'b','linewidth',2)
hold on
plot(clock, stateNon1w0psi(:,5),'r','linewidth',2)
plot(clock, stateLin1w0psi(:,4),'--r','linewidth',2)
plot(clock, stateLin1w0psi(:,5),'--b','linewidth',2)
xlabel 'Time[s]'
ylabel 'rad'
legend '\phi true' '\theta true' '\phi simulate' '\theta simulate'

f6 = figure(6);clf
f6.Position = figPos;
subplot(2,1,1)
plot(clock, stateNon1w1psi(:,1),'b','linewidth',2)
hold on
plot(clock, stateNon1w1psi(:,2),'r','linewidth',2)
plot(clock, stateLin1w1psi(:,1),'--b','linewidth',2)
plot(clock, stateLin1w1psi(:,2),'--r','linewidth',2)
hold off
xlabel 'Time[s]'
ylabel 'm'
legend 'X true' 'Y true' 'X simulate' 'Y simulate'

subplot(2,1,2)
plot(clock, stateNon1w1psi(:,4),'b','linewidth',2)
hold on
plot(clock, stateNon1w1psi(:,5),'r','linewidth',2)
plot(clock, stateLin1w1psi(:,4),'--r','linewidth',2)
plot(clock, stateLin1w1psi(:,5),'--b','linewidth',2)
xlabel 'Time[s]'
ylabel 'rad'
legend '\phi true' '\theta true' '\phi simulate' '\theta simulate'


%% We estimate the transfer function from phidot to ay

Gs = tf([0 param.g(3)*param.lambda1/param.m],[1 param.lambda1/param.m 0]);
Gd = c2d(Gs,dt,'zoh');
[bhat,ahat] = tfdata(Gd,'v');

% Simulation using A,B,C,D
A = [0 0 0 0; 0 0 0 0;0 param.g(3) -param.lambda1/param.m 0;...
    -param.g(3) 0 0 -param.lambda1/param.m];
B = [1 0 0 0;0 1 0 0]';
C = [0 0 -param.lambda1/param.m 0;0 0 0 -param.lambda1/param.m];
D = 0;
Gss = ss(A,B,C,D);


% % No wind noise and sensor fault
% y0 = acc(:,1);
% u0 = gyro(:,2);
% dat0 = estData(y0,u0,dt);
% [~,Z] = createZ(dat0,[2 2 1]);
%
% wind with no psi
yax1 = IMU1w0psi(:,1);
uq1 = IMU1w0psi(:,5);
yay1 = IMU1w0psi(:,2);
up1 = -IMU1w0psi(:,4);
datqax1 = estData(yax1,uq1);
yaxhat11 = filter(bhat,ahat,uq1);
rate1w0psi = IMU1w0psi(:,4:5);
acc1w0psi = IMU1w0psi(:,1:2);
[yaxh1,th1,xh1] = lsim(Gss,rate1w0psi,clock);

% [Y,Phi] = dat1.createZ([2 2 1]);
% th = (Z'*Phi)\(Z'*Y);
% bhat12 = [0 th(3:4)'];
% ahat12 = [1 th(1:2)'];
% yhat12 = filter(bhat12,ahat12,u1);

pm = [clock' IMU1w0psi(:,4)];
qm = [clock' IMU1w0psi(:,5)];
rm = [clock' IMU1w0psi(:,6)];

% We run the time-varying sensor to sensor model
sim('sensor2sensor.slx')
asim1w0psi = asim;


% Wind with psi change
yax2 = IMU1w1psi(:,1);
uq2 = IMU1w1psi(:,5);
yay2 = IMU1w1psi(:,2);
up2 = -IMU1w1psi(:,4);
datqax2 = estData(yax2,uq2);
yaxhat21 = filter(bhat,ahat,uq2);
rate1w1psi = IMU1w1psi(:,4:5);
acc1w1psi = IMU1w1psi(:,1:2);
[yaxh2,th2,xh2] = lsim(Gss,rate1w1psi,clock);

pm = [clock' IMU1w1psi(:,4)];
qm = [clock' IMU1w1psi(:,5)];
rm = [clock' IMU1w1psi(:,6)];

% We run the time-varying sensor to sensor model
sim('sensor2sensor.slx')
asim1w1psi = asim;



% We use the prediction error
Gs = -tf(param.g(3)*param.lambda1/param.m,[1 param.lambda1/param.m 0]);
Gd = c2d(Gs,param.dt,'zoh');
[bd,ad] = tfdata(Gd,'v');
th0 = [ad(2:end) bd(2:end)]';
% th0 = [-2+param.lambda1/param.m*param.dt 1-param.lambda1/param.m*param.dt,...
%     -param.lambda1*param.g(3)/param.m*param.dt^2]';
% th0 = th0 + 0.0*randn(3,1).*th0;
bhat = [0 th0(3:end)'];
ahat = [1 th0(1:2)'];

% INIT THE IDPOLY model for m = 0.5 kg, true mass
mARX = idpoly(ahat,bhat,[],[],[],[],param.dt);
mOE = idpoly;
mOE.f = ahat;mOE.b = bhat;mOE.Ts = param.dt;
acc0w0psi = IMU0w0psi(:,1:2);
rate0w0psi = IMU0w0psi(:,4:5);

% Data from the pitch rate to ax acc
dat1w0psi1 = iddata(acc1w0psi(:,1),rate1w0psi(:,2),param.dt);
dat1w1psi1 = iddata(acc1w1psi(:,1),rate1w1psi(:,2),param.dt);

dat0w0psi1 = iddata(acc0w0psi(:,1),rate0w0psi(:,2),param.dt);

% Data from the roll rate to ay acc
dat1w0psi2 = iddata(acc1w0psi(:,2),-rate1w0psi(:,1),param.dt);
dat1w1psi2 = iddata(acc1w1psi(:,2),-rate1w1psi(:,1),param.dt);

dat0w0psi2 = iddata(acc0w0psi(:,2),-rate0w0psi(:,1),param.dt);


OPT = peOptions('InitialCondition','z');
EqaxARX_0w0psi = pe(mARX,dat0w0psi1,1,OPT);
EqaxOE_0w0psi = pe(mOE,dat0w0psi1,1,OPT);

EqaxARX_1w0psi = pe(mARX,dat1w0psi1,1,OPT);
EqaxOE_1w0psi = pe(mOE,dat1w0psi1,1,OPT);
EqaxARX_1w1psi = pe(mARX,dat1w1psi1,1,OPT);
EqaxOE_1w1psi = pe(mOE,dat1w1psi1,1,OPT);

EpayARX_0w0psi = pe(mARX,dat0w0psi2,1,OPT);
EpayOE_0w0psi = pe(mOE,dat0w0psi2,1,OPT);

EpayARX_1w0psi = pe(mARX,dat1w0psi2,1,OPT);
EpayOE_1w0psi = pe(mOE,dat1w0psi2,1,OPT);
EpayARX_1w1psi = pe(mARX,dat1w1psi2,1,OPT);
EpayOE_1w1psi = pe(mOE,dat1w1psi2,1,OPT);

% Check using regression model for m = 1.0 kg
% th1 = [-2+param.lambda1/(2*param.m)*param.dt 1-param.lambda1/(2*param.m)*param.dt,...
%     -param.lambda1*param.g(3)/(2*param.m)*param.dt^2]';
Gs = -tf(param.g(3)*param.lambda1/(2*param.m),[1 param.lambda1/(2*param.m) 0]);
Gd = c2d(Gs,param.dt,'zoh');
[bd,ad] = tfdata(Gd,'v');
th1 = [ad(2:end) bd(2:end)]';
bhat = [0 th1(3:end)'];
ahat = [1 th1(1:2)'];

% Phi1 = [-y1(2:end-1) -y1(1:end-2) u1(2:end-1) u1(1:end-2)];
Phipay1 = [-yay1(2:end-1) -yay1(1:end-2) up1(2:end-1) up1(1:end-2)];
Phiqax1 = [-yax1(2:end-1) -yax1(1:end-2) uq1(2:end-1) uq1(1:end-2)];
Yqax1 = yax1(3:end);
Ypay1 = yay1(3:end);

eqax1w0psi = Yqax1 - Phiqax1*th1;
epay1w0psi = Ypay1 - Phipay1*th1;

% Phi2 = [-y2(2:end-1) -y2(1:end-2) u2(2:end-1) u2(1:end-2)];
Phipay2 = [-yay2(2:end-1) -yay2(1:end-2) up2(2:end-1) up2(1:end-2)];
Phiqax2 = [-yax2(2:end-1) -yax2(1:end-2) uq2(2:end-1) uq2(1:end-2)];
Yqax2 = yax2(3:end);
Ypay2 = yay2(3:end);

eqax1w1psi = Yqax2 - Phiqax2*th1;
epay1w1psi = Ypay2 - Phipay2*th1;

disp('Done residual')

% We compute the cost function based on the residual
Nw = 500;

% From no wind no psi
p0 = [zeros(Nw,1); -IMU0w0psi(:,4)];
q0 = [zeros(Nw,1); IMU0w0psi(:,5)];

ax0 = [zeros(Nw,1); IMU0w0psi(:,1)];
ay0 = [zeros(Nw,1); IMU0w0psi(:,2)];

Phi0pax1 = [-ax0(2:end-1) -ax0(1:end-2) p0(2:end-1) p0(1:end-2)];
Y0pax1 = ax0(3:end);

Phi0pax2 = [-ax0(2:end-1) -ax0(1:end-2) p0(1:end-2)];
Y0pax2 = ax0(3:end);

% From wind and no psi
p01w0psi = [zeros(Nw,1); IMU1w0psi(:,7)];
q01w0psi = [zeros(Nw,1); IMU1w0psi(:,8)];
ax01w0psi = filter(bhat,ahat,q01w0psi);
ay01w0psi = filter(bhat,ahat,p01w0psi);

% From wind and psi
p01w1psi = [zeros(Nw,1); IMU1w1psi(:,7)];
q01w1psi = [zeros(Nw,1); IMU1w1psi(:,8)];
ax01w1psi = filter(bhat,ahat,q01w1psi);
ay01w1psi = filter(bhat,ahat,p01w1psi);

dat0w0psi1 = iddata(ax0, q0, param.dt);
dat0w0psi2 = iddata(ay0, -p0, param.dt);

E01 = pe(mARX,dat0w0psi1, 1, OPT);
E02 = pe(mARX, dat0w0psi2, 1, OPT);

% CreateZ for model
Zqax0 = [-ax0(2:end-1) -ax0(1:end-2) q0(2:end-1) q0(1:end-2)];
Zpay0 = [-ay0(2:end-1) -ay0(1:end-2) p0(2:end-1) p0(1:end-2)];

Zqax1w0psi = [-ax01w0psi(2:end-1) -ax01w0psi(1:end-2) q01w0psi(2:end-1) q01w0psi(1:end-2)];
Zpay1w0psi = [-ay01w0psi(2:end-1) -ay01w0psi(1:end-2) p01w0psi(2:end-1) p01w0psi(1:end-2)];

Zqax1w1psi = [-ax01w1psi(2:end-1) -ax01w1psi(1:end-2) q01w1psi(2:end-1) q01w1psi(1:end-2)];
Zpay1w1psi = [-ay01w1psi(2:end-1) -ay01w1psi(1:end-2) p01w1psi(2:end-1) p01w1psi(1:end-2)];

datcost = estData(zeros(N,1),EqaxARX_1w0psi.y);
% [~,Zcost] = datcost.createZ([0 Nw 0]);
YqaxARX1 = [zeros(Nw,1); EqaxARX_1w0psi.y(3:end)];
YqaxARX2 = [zeros(Nw,1); EqaxARX_1w1psi.y(3:end)];

Y2mqaxARX1 = [zeros(Nw,1); eqax1w0psi];
Y2mqaxARX2 = [zeros(Nw,1); eqax1w1psi];

YqaxOE1 = [zeros(Nw,1); EqaxOE_1w0psi.y(3:end)];
YqaxOE2 = [zeros(Nw,1); EqaxOE_1w1psi.y(3:end)];

YpayARX1 = [zeros(Nw,1); EpayARX_1w0psi.y(3:end)];
YpayARX2 = [zeros(Nw,1); EpayARX_1w1psi.y(3:end)];

Y2mpayARX1 = [zeros(Nw,1); epay1w0psi];
Y2mpayARX2 = [zeros(Nw,1); epay1w1psi];

YpayOE1 = [zeros(Nw,1); EpayOE_1w0psi.y(3:end)];
YpayOE2 = [zeros(Nw,1); EpayOE_1w1psi.y(3:end)];

Vu = pdotref(:,1);
Vv = pdotref(:,2);

Vuu = [zeros(Nw,1); Vu];

Lambda = ones(Nw+2,1);
for kk_p = 1:length(Lambda)
    Lambda(kk_p) = 0.999^(length(Lambda)-kk_p+1);
end
Lambda = sqrt(Lambda);

for ii = 1:size(Zqax0,1)-Nw-1
    ZZqax0 = Zqax0(ii:Nw+ii+1,1:end);
    ZZpay0 = Zpay0(ii:Nw+ii+1,1:end);
    
    ZZqax1w0psi = Zqax1w0psi(ii:Nw+ii+1,1:end);
    ZZpay1w0psi = Zpay1w0psi(ii:Nw+ii+1,1:end);
    
    ZZqax1w1psi = Zqax1w1psi(ii:Nw+ii+1,1:end);
    ZZpay1w1psi = Zpay1w1psi(ii:Nw+ii+1,1:end);
    
    
    JqaxARX1(:,ii) = norm(ZZqax0.*(YqaxARX1(ii:Nw+ii+1,1).*Lambda));
    JqaxARX2(:,ii) = norm(ZZqax0.*YqaxARX2(ii:Nw+ii+1,1));
    
    JpayARX1(:,ii) = norm(ZZpay0.*YpayARX1(ii:Nw+ii+1,1));
    JpayARX2(:,ii) = norm(ZZpay0.*YpayARX2(ii:Nw+ii+1,1));
    
    JqaxOE1(:,ii) = mean(YqaxOE1(ii:Nw+ii,1).^2);
    JqaxOE2(:,ii) = mean(YqaxOE2(ii:Nw+ii,1).^2);
    
    % ZpayARX1(:,ii) = mean(YpayARX1(ii:Nw+ii,1).^2);
    % ZpayARX2(:,ii) = mean(YpayARX2(ii:Nw+ii,1).^2);
    
    JpayOE1(:,ii) = mean(YpayOE1(ii:Nw+ii,1).^2);
    JpayOE2(:,ii) = mean(YpayOE2(ii:Nw+ii,1).^2);
    
    Jqax1w0psi(ii,:) = norm(ZZqax0'*(YqaxARX1(ii:Nw+ii+1,1).*Lambda));
    Jqax1w1psi(ii,:) = norm(ZZqax0'*(YqaxARX2(ii:Nw+ii+1,1).*Lambda));
    
    Jqax1w0psi1(ii,:) = norm(ZZqax1w0psi'*(YqaxARX1(ii:Nw+ii+1,1).*Lambda));
    Jqax1w1psi1(ii,:) = norm(ZZqax1w1psi'*(YqaxARX2(ii:Nw+ii+1,1).*Lambda));
    
    
    Jpay1w0psi(ii,:) = norm(ZZpay0'*(YpayARX1(ii:Nw+ii+1,1).*Lambda));
    Jpay1w1psi(ii,:) = norm(ZZpay0'*(YpayARX2(ii:Nw+ii+1,1).*Lambda));
    
    Jpay1w0psi1(ii,:) = norm(ZZpay1w0psi'*(YpayARX1(ii:Nw+ii+1,1).*Lambda));
    Jpay1w1psi1(ii,:) = norm(ZZpay1w1psi'*(YpayARX2(ii:Nw+ii+1,1).*Lambda));
    
    J2mqax1w0psi(ii,:) = norm(ZZqax0'*(Y2mqaxARX1(ii:Nw+ii+1,1).*Lambda));
    J2mqax1w1psi(ii,:) = norm(ZZqax0'*(Y2mqaxARX2(ii:Nw+ii+1,1).*Lambda));
    J2mpay1w0psi(ii,:) = norm(ZZpay0'*(Y2mpayARX1(ii:Nw+ii+1,1).*Lambda));
    J2mpay1w1psi(ii,:) = norm(ZZpay0'*(Y2mpayARX2(ii:Nw+ii+1,1).*Lambda));
    
end


disp('Done cost function')

% [Y,Phi] = dat2.createZ([2 2 1]);
% th = (Z'*Phi)\(Z'*Y);
% bhat22 = [0 th(3:4)'];
% ahat22 = [1 th(1:2)'];
% yhat22 = filter(bhat22,ahat22,u1);

% % No wind and sensor fault
% y2 = acc0(:,1);
% u2 = gyro0(:,2);
% yy3 = acc0(:,1:2) + param.lambda1/param.m*windm(:,4:5);
% uu3 = gyro(:,1:2);
% [yh3,th3,xh3] = lsim(Gss,uu3,clock);

% PLOT THE FIGURE
close all
f7 = figure(7);clf
f7.Position = [ -1900         680         560         420];
h1 = subplot(2,1,1);
h1.FontSize = 11;
hold on
% plot(clock,yy1(:,1)-yh1(:,1));
% plot(clock,yy2(:,1)-yh2(:,1),'--')
plot(clock,EqaxARX_1w0psi.y,'-','linewidth',1.5)
plot(clock,EqaxOE_1w0psi.y,'-.','linewidth',1.5)
plot(clock,asim1w0psi(:,1)+acc1w0psi(:,1),'--')

plot(clock,EqaxARX_1w1psi.y,'-','linewidth',1.5)
plot(clock,EqaxOE_1w1psi.y,'-.','linewidth',1.5)
plot(clock,asim1w1psi(:,1)+acc1w1psi(:,1))
hold off
xlabel('Time[s]','FontSize',18)
ylabel({'Residual $r_1$'},'Interpreter','latex','FontSize',18)
legend('without \psi change + ARX model, ','without \psi change + OE model','without \psi change + varying',...
    'with \psi change + ARX model','with \psi change + OE model',...
    'with \psi change + varying','FontSize',16)
% ylim([-1 2])

h2 = subplot(2,1,2);
h2.FontSize = 11;
hold on
% plot(clock,yy1(:,2)-yh1(:,2));
% plot(clock,yy2(:,2)-yh2(:,2),'--')
plot(clock,EpayARX_1w0psi.y,'-','linewidth',1.5)
plot(clock,EpayOE_1w0psi.y,'-.','linewidth',1.5)
plot(clock,asim1w0psi(:,2)+acc1w0psi(:,2),'--')

plot(clock,EpayARX_1w1psi.y,'-','linewidth',1.5)
plot(clock,EpayOE_1w1psi.y,'-.','linewidth',1.5)
plot(clock,asim1w1psi(:,2)+acc1w1psi(:,2))
hold off
xlabel('Time[s]','FontSize',18)
ylabel({'Residual $r_2$'},'Interpreter','latex','FontSize',18)
legend('without \psi change + ARX model, ','without \psi change + OE model','without \psi change + varying',...
    'with \psi change + ARX model','with \psi change + OE model',...
    'with \psi change + varying','FontSize',16)
% ylim([-2 2])

f8 = figure(8);clf
f8.Position = [-1250         600         560         420];
h1 = subplot(2,1,1);
h1.FontSize = 11;
hold on
% plot(clock,yy1(:,1)-yh1(:,1));
% plot(clock,yy2(:,1)-yh2(:,1),'--')
plot(clock(4:end),Jqax1w0psi,'-','linewidth',1.5)
plot(clock(4:end),Jqax1w1psi,'-.','linewidth',1.5)
% plot([40 40],[0 1],'-','linewidth',1.5)
hold off
xlabel('Time[s]','FontSize',18)
ylabel({'The IV cost function $J_{IV1}$'},'Interpreter','latex','FontSize',18)
legend('without \psi change','with \psi change ','FontSize',16)
% ylim([-.1 1])
title('Model: pitch rate to x-axis acceleration: mass change at 40s','FontSize',18)
% ylim([0 0.2])

h2 = subplot(2,1,2);
h2.FontSize = 11;
hold on
plot(clock(4:end),Jpay1w0psi,'-','linewidth',1.5)
plot(clock(4:end),Jpay1w1psi,'-.','linewidth',1.5)
% plot(clock(4:end),sign(mt(4:end)-0.5),'-','linewidth',2)
hold off
xlabel('Time[s]','FontSize',18)
ylabel({'The IV cost function $J_{IV2}$'},'Interpreter','latex','FontSize',18)
legend('without \psi change, ','with \psi change','FontSize',16)
title('Model: roll rate to y-axis acceleration: mass change at 40s','FontSize',18)
% ylim([-.1 1])

f9 = figure(9);clf
f9.Position = [-650   600   560   420];
h1 = subplot(2,1,1);
h1.FontSize = 11;
hold on
% plot(clock,yy1(:,1)-yh1(:,1));
% plot(clock,yy2(:,1)-yh2(:,1),'--')
plot(clock(4:end),Jqax1w0psi,'-','linewidth',1.5)
plot(clock(4:end),J2mqax1w0psi,'-.','linewidth',1.5)
% plot([40 40],[0 1],'-','linewidth',1.5)
hold off
xlabel('Time[s]','FontSize',18)
ylabel({'The IV cost function $J_{IV1}$'},'Interpreter','latex','FontSize',18)
legend('with m=0.5 kg','with m=1 kg','FontSize',16)
% ylim([-.1 1])
title('Model: pitch rate to x-axis acceleration: without \psi change','FontSize',18)
% ylim([0 0.2])

h2 = subplot(2,1,2);
h2.FontSize = 11;
hold on
plot(clock(4:end),Jqax1w0psi - J2mqax1w0psi,'-','linewidth',1.5)
% plot(clock(4:end),sign(mt(4:end)-0.5),'-','linewidth',2)
plot([clock(4) clock(end)],[0 0],'linewidth',2)
% hold off
xlabel('Time[s]','FontSize',18)
ylabel({'Difference $J_{IV1}$'},'Interpreter','latex','FontSize',18)
% legend('without \psi change, ','with \psi change','FontSize',16)
title('Difference of J_{IV1} between m=0.5kg and m=1kg','FontSize',18)
hh=get(gca,'Children'); % grab all the axes handles at once
legendstr={'','Zero threshold'};
legend(hh(1),legendstr(2), 'Fontsize',18,'Location', 'southeast', 'Orientation','horizontal') % see the second label is skipped withoud changing shape
hold off
% legend '' 'Zero threshold'
% ylim([-5e-5 1e-5])

f10 = figure(10);clf
f10.Position = [ -1850          90         560         420];
h1 = subplot(2,1,1);
h1.FontSize = 11;
hold on
% plot(clock,yy1(:,1)-yh1(:,1));
% plot(clock,yy2(:,1)-yh2(:,1),'--')
plot(clock(4:end),Jqax1w1psi,'-','linewidth',1.5)
plot(clock(4:end),J2mqax1w1psi,'-.','linewidth',1.5)
% plot([40 40],[0 1],'-','linewidth',1.5)
hold off
xlabel('Time[s]','FontSize',18)
ylabel({'The IV cost function $J_{IV1}$'},'Interpreter','latex','FontSize',18)
legend('with m=0.5 kg','with m=1 kg','FontSize',16)
% ylim([-.1 1])
title('Model: pitch rate to x-axis acceleration: with \psi change','FontSize',18)
% ylim([0 0.2])

h2 = subplot(2,1,2);
h2.FontSize = 11;
hold on
plot(clock(4:end),Jqax1w1psi - J2mqax1w1psi,'-','linewidth',1.5)
% plot(clock(4:end),sign(mt(4:end)-0.5),'-','linewidth',2)
plot([clock(4) clock(end)],[0 0],'linewidth',2)
% hold off
xlabel('Time[s]','FontSize',18)
ylabel({'Difference $J_{IV1}$'},'Interpreter','latex','FontSize',18)
% legend('without \psi change, ','with \psi change','FontSize',16)
title('Difference of J_{IV1} between m=0.5kg and m=1kg','FontSize',18)
hh=get(gca,'Children'); % grab all the axes handles at once
legendstr={'','Zero threshold'};
legend(hh(1),legendstr(2), 'Fontsize',18,'Location', 'southeast', 'Orientation','horizontal') % see the second label is skipped withoud changing shape
hold off
% ylim([-5e-5 1e-5])


f11 = figure(11);clf;
f11.Position = [ -1250    90   560   420];
h1 = subplot(2,1,1);
h1.FontSize = 11;
hold on
plot(clock(4:end),Jqax1w0psi,'-','linewidth',1.5)
% hold off
xlabel('Time[s]','FontSize',18)
ylabel({'$J_{IV1}(m=0.5)$'},'Interpreter','latex','FontSize',18)
% legend('with m=0.5 kg','FontSize',16)
% ylim([-.1 1])
title('Model: pitch rate to x-axis acceleration: without \psi change','FontSize',18)
% ylim([0 0.2])

h2 = subplot(2,1,2);
h2.FontSize = 11;
hold on
plot(clock(4:end),Jqax1w1psi,'-','linewidth',1.5)
% hold off
xlabel('Time[s]','FontSize',18)
ylabel({'$J_{IV1}(m=0.5)$'},'Interpreter','latex','FontSize',18)
% legend('without \psi change, ','with \psi change','FontSize',16)
title('Model: pitch rate to x-axis acceleration: with \psi change','FontSize',18)
% legend('with m=0.5 kg','FontSize',16)
hold off
% ylim([-5e-5 1e-5])

%%
TT = 30/param.dt:50/param.dt;
figure(9);clf
plot(clock(30/param.dt+4:50/param.dt+4),Rqax1w0psi(30/param.dt:50/param.dt),'-','linewidth',1.5)
hold on
plot(clock(30/param.dt+4:50/param.dt+4),J2mqax1w0psi(30/param.dt:50/param.dt),'-.','linewidth',1.5)
xlim([30 50]);
ylim([0 0.2])
set(gca,'YTick',[0 0.1 0.2],'XTick',[30 35 40 45 50])
set(gca,'FontSize',18)

figure(10);clf
plot(clock(30/param.dt+4:50/param.dt+4),Rpay1w0psi(30/param.dt:50/param.dt),'-','linewidth',1.5)
hold on
plot(clock(30/param.dt+4:50/param.dt+4),J2mpay1w0psi(30/param.dt:50/param.dt),'-.','linewidth',1.5)
xlim([30 50]);
ylim([0 0.2])
set(gca,'YTick',[0 0.1 0.2],'XTick',[30 35 40 45 50])
set(gca,'FontSize',18)

%%
figure(7);clf
h1 = subplot(2,1,1);
h1.FontSize = 11;
hold on;
plot(clock, uq1,clock,u2,'linewidth',2)
xlabel 'Time[s]'
ylabel({'Pitch rate $\dot{\phi}$'},'Interpreter','latex')
legend('without wind','with wind')

h2 = subplot(2,1,2);
h2.FontSize = 11;
hold on
plot(clock,y1,clock,y2,'linewidth',2)
xlabel 'Time[s]'
ylabel('Longitudinal acc a_y')
legend('without wind','with wind')


figure(8);clf
h1 = subplot(2,1,1);
h1.FontSize = 11;
hold on
plot(clock,acc1w0psi(:,1)-yh1(:,1));
plot(clock,yy3(:,1)-yh3(:,1),'--')
hold off
xlabel 'Time[s]'
ylabel({'Residual $r_1$'},'Interpreter','latex')
legend('without wind','with acc bias')

h2 = subplot(2,1,2);
h2.FontSize = 11;
hold on
plot(clock,acc1w0psi(:,2)-yh1(:,2));hold on;
plot(clock,yy3(:,2)-yh3(:,2),'--')
hold off
xlabel 'Time[s]'
ylabel({'Residual $r_2$'},'Interpreter','latex')
legend('without wind','with acc bias')

figure(9);clf
h1 = subplot(2,1,1);
h1.FontSize = 11;
hold on;
plot(clock, uq1,clock,u2,'linewidth',2)
xlabel 'Time[s]'
ylabel({'Pitch rate $\dot{\phi}$'},'Interpreter','latex')
legend('without wind','with wind')

h2 = subplot(2,1,2);
h2.FontSize = 11;
hold on
plot(clock,y1,clock,y2,'linewidth',2)
xlabel 'Time[s]'
ylabel('Longitudinal acc a_y')
legend('without wind','with wind')



