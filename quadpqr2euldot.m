function [Rpqr,euldot] = quadpqr2euldot(euler,pqr)
%quadPqr2euldot from pqr to euler dot
%   If we want to call the function from SIMULINK, the output should be a
%   vector
Rpqr = quadRpqr2euldot(euler);
euldot = Rpqr*pqr;
Rpqr = reshape(Rpqr,9,1);
end

