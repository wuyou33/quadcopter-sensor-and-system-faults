function S = quadSskew(omega)
%quadSskew Function to compute the skew symmetric of angular velocities
%   Detailed explanation goes here

if length(omega) == 1
    S = [0 -omega;omega 0];
elseif length(omega) == 3
    S = [0 -omega(3) omega(2);...
        omega(3) 0 -omega(1);...
        -omega(2) omega(1) 0];
end
