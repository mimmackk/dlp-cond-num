function k = cond_num(P, eigval, x, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% Calculates the condition number for matrix polynomial for a given
% eigenvalue and left and right eigenvectors.
%
% INPUT
% P: cell array of matrix coefficients
% eigval: eigenvalue to evaluate at
% x: right eigenvector
% y: left eigenvector
%
% OUTPUT
% k: condition number of P at given eigenvalue, eigenvectors
%
% AUTHORS
% Written by ____
% Minor modifications by Kayden Mimmack & Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets digit precision higher for more accurate computations
digits(40);

k = length(P) - 1;

syms d;
kappa1 = symfun(0, d); 
kappa2 = symfun(derivP(P), d);

for j = 0 : k
    kappa1 = kappa1 + abs(d)^j * norm(P{j + 1});
end

% Construct numerator and denominator separately, then divide
kappa1 = kappa1 * norm(y) * norm(x);
kappa2 = abs(d) * abs(ctranspose(y) * kappa2 * x);
k = kappa1 / kappa2;

% Evaluate at eigenvalue
k = k(eigval);