function [xref, yref, zref] = quadtrapezoidal(clock)
% quadtrapezoidal Path generator using trapezoidal, bandlimited speed
phimax = 0.1;
lambda1 = 0.36;
m = 0.5;
g = 9.81;
vmax = m*g*phimax/lambda1;
vzmax = 2;

thover = 10;
tx1 = 25;
ty1 = 40;
tx2 = 55;
ty2 = 70;
T = 5;
amax = vmax/T;

if clock<=thover && clock>=0 % From ground to hovering position
    
    yref = 0;xref = 0;
    zref = vzmax*clock;
    if zref>=0.8
        zref = 0.8;
    end
end
if clock>=thover && clock <=tx1
    zref = 0.8;
    if clock>=thover && clock <=thover+T
        xref = 0.5*amax*(clock-thover)^2;
        yref = 0;
    end
    if clock>=thover+T && clock <=tx1-T
        xref = 0.5*amax*T^2 + vmax*(clock-(thover+T));
        yref = 0;
    end
    if clock>=tx1-T && clock <=tx1
        xref = 0.5*amax*T^2 + vmax*(tx1-thover-2*T)+...
            vmax*(clock-tx1+T)-0.5*amax*(clock-tx1+T)^2;
        yref = 0;
    end
end
if clock>=tx1
    xref = vmax*(tx1-thover-T);
    yref = 0;
    zref = 0.8;
end
end