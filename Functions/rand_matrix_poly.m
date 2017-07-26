function P = rand_matrix_poly(k, n, r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% rand_matrix_poly generates a random matrix polynomial of degree k, size n
%
% INPUT
% k: desired degree
% n: desired size (square)
% r: range of values desired in P: entries between -r and r
%
% OUTPUT
% P: Symbolic matrix poly (cell array of matrix coefficents A_0 to A_k)
%
% AUTHORS
% Written by ____
% Modifications by Kayden Mimmack & Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = cell(1, k + 1);

for i = 1 : (k + 1)
    
    % Generate an n x n random matrix with entries between -r and r
    P{i} = sym(-r + 2 * r * rand(n), 'f');
    
end
    