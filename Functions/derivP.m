function dP = derivP(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% derivP calcualtes the symbolic derivative of a matrix polynomial
%
% INPUT
% P: Symbolic matrix poly (cell array of matrix coefficents A_0 to A_k)
%
% OUTPUT
% dP: Symbolic derivative of P
%
% AUTHORS
% Written by Lucas Medina
% Modifications by Kayden Mimmack & Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = length(P) - 1;

syms d;
f = symfun(1, d) * P{1};

for i = 1 : k
    f = f + symfun(d^i, d) * P{i + 1};
end

dP = diff(f);
        