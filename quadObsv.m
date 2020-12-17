function [O,rankO] = quadObsv(f,varargin)
%quadObsv Compute the observability matrix for time-varying yaw models
if nargin == 4 % one h
    h1 = varargin{1};
    h2 = [];
    psivec = varargin{2};
    rvec = varargin{3};
elseif nargin == 5 % two h
    h1 = varargin{1};
    h2 = varargin{2};
    psivec = varargin{3};
    rvec = varargin{4};
end
nr = length(rvec);
np = length(psivec);
if nr ~=np
    error('wrong length')
end
% Compute the observability matrix
if isempty(h2)
    C1 = h1(psivec(1),rvec(1));
else 
    C1 = [h1(psivec(1),rvec(1));h2(psivec(1),rvec(1))];
end
A1 = f(psivec(1),rvec(1));
[nc, mc] = size(C1);
O = zeros(nc*nr,mc);
Cii = zeros(nc,mc,nr);
Aii = zeros(mc,mc,nr);
for ii = 1:nr
    if isempty(h2)
        hii = h1(psivec(ii),rvec(ii));
    else 
        hii = [h1(psivec(ii),rvec(ii));h2(psivec(ii),rvec(ii))];
    end
    Cii(:,:,ii) = hii;
    Aii(:,:,ii) = f(psivec(ii),rvec(ii));
    if ii == 1
        O(1:nc,:) = C1;
    else
        Ajj = eye(mc);
        for jj = 1:ii-1
            Ajj = Ajj*squeeze(Aii(:,:,ii-jj));
        end
        O((ii-1)*nc+1:ii*nc,:) = squeeze(Cii(:,:,ii))*Ajj;
        
    end
end
rankO = rank(O'*O);

