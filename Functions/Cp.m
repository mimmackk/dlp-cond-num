function [Cp, iAk] = Cp(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% Cp generates linearization of P C(P)
%
% INPUT
% P: Cell array of matrix coefficients
%
% OUTPUT
% Cp: C(P) linearization associated with P
% iAk: inverse of A_k
%
% AUTHORS
% Written by ____
% Minor modifications by Kayden Mimmack & Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ak = P{length(P)};
n = size(Ak);
k = length(P) - 1;

% Generate inverse of A_k
iAk = Ak \ eye(n);

% Initialize Cp to be symbolic matrix of size nk full of zeros
Cp = zeros(n * k);
Cp = sym(Cp, 'f');

% Place -inv(Ak) * A_{i} for i = k-1, k-2, ..., 0 on left partitions
for j = 0 : (k - 1) 
    Cp((1 + j * n) : (j * n + n) , 1 : n) = -iAk * P{length(P) - (j + 1)};
end

j = 1;
z = n + 1;

% Put identity matrices in the top right partition
while z + 1 <= (n * k)
    Cp(j : (j + (n - 1)), z : (z + (n - 1))) = eye(n);
    j = j + n;
    z = z + n;
end