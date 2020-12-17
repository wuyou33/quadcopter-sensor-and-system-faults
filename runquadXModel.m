%%
clc
clear all
close all

dt = 1/200;
tfinal = 100;

t = 0:dt:tfinal;t = t';
g = 9.81;
lambda = 0.36;
T = dt;
Gs = tf(g*lambda/0.5,[1 lambda/0.5 0]);
[bd1, ad1] = tfdata(c2d(Gs,dt,'zoh'),'v');
[bd2, ad2] = tfdata(c2d(Gs,dt,'foh'),'v');
[bd3, ad3] = tfdata(c2d(Gs,dt,'tustin'),'v');


vfinal = 0;tmidd = 0;
TxL1 = [0 tmidd]; TxL2 = tmidd:dt:tfinal;
VxL1 = [0 vfinal]; VxL2 = vfinal + 0.5*sin(2*pi/2.5*(TxL2-tmidd)).*cos(2*pi/2*(TxL2-tmidd)) + ...
    0.5*sin(2*pi/2*(TxL2-tmidd)).*cos(2*pi/1.5*(TxL2-tmidd)) + ...
    0.5*sin(2*pi/1.5*(TxL2-tmidd)).*cos(2*pi/1*(TxL2-tmidd));

ref = [[TxL1 TxL2]' [VxL1 VxL2]'];
% ref = [t filter(0.5,[1 -0.5],randn(size(t)))];
% Create Instrument

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
Hud = c2d(Hu,dt);
[bHu,aHu] = tfdata(Hud,'v');


for ii = 0:2

    switch(ii)
        case 0
            tm = [0 50 52 tfinal];
            m = [0.5 0.5 0.5 0.5];
            mass = [tm' m'];
            
            
            ea = [t sqrt(0.)*randn(size(t))];
            ep = [t sqrt(0.)*randn(size(t))];
            windnoise = [t 0+0*filter(0.05,[1 -0.95],sqrt(1)*randn(size(t)))];
            sim('quadXmodel.slx')
            
            emodel00 = emodel;
            ax00 = axm;
            p00 = pm;
            disp('Zero case')

        case 1
            tm = [0 50 52 tfinal];
            m = [0.5 0.5 1 1];
            mass = [tm' m'];
            
            
            ea = [t sqrt(0.)*randn(size(t))];
            ep = [t sqrt(0.)*randn(size(t))];
            windnoise = [t 1+filter(bHu,aHu,sqrt(10)*randn(size(t)))];
            sim('quadXmodel.slx')
            
            emodel0 = emodel;
            ax0 = axm;
            p0 = pm;
            disp('First case')

        case 2
            tm = [0 50 52 tfinal];
            m = [0.5 0.5 1 1];
            mass = [tm' m'];
            
            ea = [t sqrt(1e-6)*randn(size(t))];
            ep = [t sqrt(1e-12)*randn(size(t))];
            windnoise = [t 1+filter(bHu,aHu,sqrt(10)*randn(size(t)))];
            sim('quadXmodel.slx')
            
            axm = axm;
            pm = pm;
            disp('Second case')

    end
    
end
%%
g = 9.81;
T = dt;
lambda = 0.36;
Gs = tf(g*lambda/0.5,[1 lambda/0.5 0]);
[bd1, ad1] = tfdata(c2d(Gs,dt,'zoh'),'v');
[bd2, ad2] = tfdata(c2d(Gs,dt,'foh'),'v');
[bd3, ad3] = tfdata(c2d(Gs,dt,'tustin'),'v');

th1 = [ad1(2:end) bd1(2)+bd1(3)]';

fth  =  @(x,th,Pth) (th - [-2+lambda/x*T, 1-lambda/x*T, g*lambda/x*T^2]')'*inv(Pth)*(th - [-2+lambda/x*T, 1-lambda/x*T, g*lambda/x*T^2]');

%


Gs = tf(g*lambda/1,[1 lambda/1 0]);
Gd = c2d(Gs,dt,'zoh');
[bd2, ad2] = tfdata(Gd,'v');
th2 = [ad2(2:end) bd2(2)+bd2(3)]';

% ax00 = ax0;
% p00 = p0;

% ax00 = IMU0w0psi(:,1);
% p00 = -IMU0w0psi(:,5);
% 
% ax0 = IMU1w0psi(:,1);
% p0 = -IMU1w0psi(:,5);
% 
% axm = IMU1w1psi(:,1);
% pm = -IMU1w1psi(:,5);


Phi00 = [-ax00(2:end-1) -ax00(1:end-2) p00(1:end-2)];
Y00 = ax00(3:end);

Phi0 = [-ax0(2:end-1) -ax0(1:end-2) p0(1:end-2)];
Y0 = ax0(3:end);

Phim = [-axm(2:end-1) -axm(1:end-2) pm(1:end-2)];
Ym = axm(3:end);

N = length(Y0);
Nw = 500;
for kk = 1:Nw+1
lambdaForget(kk,:) = 0.999^(Nw-kk+1);
end

% Apply the Kalman filter
fkf = @(x,u) [x(1)+dt*u; x(2)+g*dt*x(1)-lambda*dt/x(3)*x(2); x(3)];
hkf = @(x,u) lambda/x(3)*x(2);

dfkf = @(x,u) [1 0 0;g*dt 1-lambda*dt/x(3) lambda*dt*x(2)*x(3)^2;0 0 1];
dhkf = @(x,u) [0 lambda/x(3) -lambda*x(2)/x(3)^2];
R = 1e-5;
Q = diag([1e-5 1e-5 1e-5]);
x = [0 0 0.25]';
P = 1000*eye(3);

thr = [-1.9964 0.9964 0.0017]';
Pr = 1000*eye(3);
for kk = 1:length(Y00)
    
    ZZ0 = Phi00(kk,:)';
    Phi02 = Phim(kk,:)';
    yy0 = Ym(kk,:);
    Pr = Pr + diag([0.01 0.01 0.001]);
    K = Pr*ZZ0*inv(1 + Phi02'*Pr*ZZ0);
    Pr = Pr - K*Phi02'*Pr;
    thr = thr + K*(yy0-Phi02'*thr);
    
    thrmc(kk,:) = thr;
    [x,P]=ekf_cao(fkf,hkf,dfkf,dhkf,x,P,axm(kk),p0(kk),Q,R);
    xmc(kk,:) = x;
end

for ii = 1:N-Nw
    Z0 = Phi00(ii:ii+Nw,:);
    YY = Y00(ii:ii+Nw,:);
    
    Phi1 = Phi0(ii:ii+Nw,:);
    Y1 = Y0(ii:ii+Nw,:);
    
    Phi2 = Phim(ii:ii+Nw,:);
    Y2 = Ym(ii:ii+Nw,:);
    
    th0 = (Z0'*Phi1)\(Z0'*Y1);
    th01 = Phi1\Y1;
    ee0 = Y1-Phi1*th0;
    Pth = inv(Z0'*Phi1)*(Z0'*ee0*ee0'*Z0)*inv(Z0'*Phi1)';    
    condmc(ii,:) = cond(Z0'*Phi1);
    th0mc(ii,:) = th0;
    
    th = (Z0'*Phi2)\(Z0'*Y2);
    thmc(ii,:) = th;
    

    
%     m = 0.5;
%     Xm = fminsearch(@(x) fth(x,th0,eye(3)), m);
%     Xmc(ii,:) = Xm;
%     
    
    sig11(ii,:) = norm(Z0'*((Y1-Phi1*th1).*lambdaForget));
    sig21(ii,:) = norm(Z0'*((Y1-Phi1*th2).*lambdaForget));
    
    sig12(ii,:) = norm(Z0'*((Y2-Phi2*th1).*lambdaForget));
    sig22(ii,:) = norm(Z0'*((Y2-Phi2*th2).*lambdaForget));

end
% tSimt = 0:dt:dt*(length(Xmc)-1);
figure(11);clf
plot(sig11)
hold on
plot(sig21)
plot(sig12,'--')
plot(sig22,'-.')
% ylim([-1e-5 1e-2])
hold off
legend 'No noise 0.5' 'No noise 1' 'Noise 0.5' 'Noise 1'

figure(12)
plot(xmc(:,3))
hold on
plot(mSim)
% plot(Xmc)
hold off






