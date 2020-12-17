function [Oout,productAt] = quadObsvRecurive(At,Ct,Opre,productAtpre)
 % quadObsvRecurive Compute the observability matrix recursively
    Oout = [Opre; Ct*productAtpre];
    productAt = At*productAtpre;
end