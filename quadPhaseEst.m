%%
clc
tt = 0:0.001:10;
phase = pi/2-rand*pi;

for mc = 1:100
    R = 1e-2;
    Q = diag([1e-6 1e-6]);
exx = sqrt(R)*randn(size(tt));
xx = sin(2*pi*0.1*tt+phase)+exx;

% Now use KF to estimate
xp = [0;0];
Pp = 1000*eye(2);
xpn = xp;
Ppn = Pp;
for ii = 1:length(tt)
    xp = eye(2)*xp;
    Pp = Pp + Q;
    
    Xpn = eye(2)*xpn;
    Ppn = Ppn+ Q;
    
    % Correction step
    Hp = [sin(2*pi*0.1*tt(ii)) cos(2*pi*0.1*tt(ii))];
    Sp = Hp*Pp*Hp'+R;
    Kp = Pp*Hp'/Sp;
    xp = xp + Kp*(xx(ii)-Hp*xp);
    Pp = Pp-Kp*Hp*Pp;
    xpm(ii,:,mc) = xp;
    stdp(ii,:,mc) = diag(Pp);
    
    
    % Correction step
    eii = sqrt(1e-3)*randn;
    Hpn = [sin(2*pi*0.1*tt(ii)+eii) cos(2*pi*0.1*tt(ii)+eii)];
    Spn = Hpn*Ppn*Hpn'+R;
    Kpn = Ppn*Hpn'/Spn;
    xpn = xpn + Kpn*(xx(ii)-Hpn*xpn);
    Ppn = Ppn-Kpn*Hpn*Ppn;
    xpnm(ii,:,mc) = xpn;
    stdpn(ii,:,mc) = diag(Ppn);
end


mc
end

%%
x1m = squeeze(xpm(:,1,:));
x2m = squeeze(xpm(:,2,:));
stdp1 = squeeze(stdp(:,1,:));
stdp2 = squeeze(stdp(:,2,:));
mm = mean(xpm,3);

m_x1m = mean(x1m');
std_x1m = std(x1m');
std_p1 = mean(stdp1');

m_x2m = mean(x2m');
std_x2m = std(x2m');
std_p2 = mean(stdp2');

tt2 = [tt tt(end:-1:1)];
std_x11 = [m_x1m-std_x1m m_x1m(end:-1:1)+std_x1m(end:-1:1)];
std_x21 = [m_x2m-std_x2m m_x2m(end:-1:1)+std_x2m(end:-1:1)];

std_x12 = [m_x1m-std_p1 m_x1m(end:-1:1)+std_p1(end:-1:1)];
std_x22 = [m_x2m-std_p2 m_x2m(end:-1:1)+std_p2(end:-1:1)];



x1mn = squeeze(xpnm(:,1,:));
x2mn = squeeze(xpnm(:,2,:));
stdp1n = squeeze(stdpn(:,1,:));
stdp2n = squeeze(stdpn(:,2,:));


m_x1mn = mean(x1mn');
std_x1mn = std(x1mn');
std_p1n = mean(stdp1n');

m_x2mn = mean(x2mn');
std_x2mn = std(x2mn');
std_p2n = mean(stdp2n');

tt2 = [tt tt(end:-1:1)];
std_x11n = [m_x1mn-std_x1mn m_x1mn(end:-1:1)+std_x1mn(end:-1:1)];
std_x21n = [m_x2mn-std_x2mn m_x2mn(end:-1:1)+std_x2mn(end:-1:1)];

std_x12n = [m_x1mn-std_p1n m_x1mn(end:-1:1)+std_p1n(end:-1:1)];
std_x22n = [m_x2mn-std_p2n m_x2mn(end:-1:1)+std_p2n(end:-1:1)];


figure(1);clf
subplot(2,1,1)
plot([0 tt(end)],cos(phase)*[1 1],'linewidth',2)
hold on
plot(tt,m_x1m,'--r','linewidth',2)
patch('Faces',(1:length(tt2)),'Vertices',[tt2' std_x11'],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(tt2)),'Vertices',[tt2' std_x12'],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.4);


plot(tt,m_x1mn,'-green','linewidth',2)
patch('Faces',(1:length(tt2)),'Vertices',[tt2' std_x11n'],'EdgeColor','none','FaceColor','black','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(tt2)),'Vertices',[tt2' std_x12n'],'EdgeColor','none','FaceColor','black','LineWidth',2,'FaceAlpha',0.4);
ylim([-1 1]+cos(phase))


subplot(2,1,2)
plot([0 tt(end)],sin(phase)*[1 1],'linewidth',2)
hold on
plot(tt,m_x2m,'--r','linewidth',2)
patch('Faces',(1:length(tt2)),'Vertices',[tt2' std_x21'],'EdgeColor','none','FaceColor','red','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(tt2)),'Vertices',[tt2' std_x22'],'EdgeColor','none','FaceColor','green','LineWidth',2,'FaceAlpha',0.4);
plot(tt,m_x2mn,'-green','linewidth',2)
patch('Faces',(1:length(tt2)),'Vertices',[tt2' std_x21n'],'EdgeColor','none','FaceColor','black','LineWidth',2,'FaceAlpha',0.3);
patch('Faces',(1:length(tt2)),'Vertices',[tt2' std_x22n'],'EdgeColor','none','FaceColor','black','LineWidth',2,'FaceAlpha',0.4);
ylim([-1 1]+sin(phase))













