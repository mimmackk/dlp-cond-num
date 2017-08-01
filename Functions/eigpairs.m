function [E, V] = eigpairs(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% Takes a matrix polynomial and returns its eigenvalues. Requires use of Cp
% linearization in order to do this.
%
% INPUT
% P: Matrix polynomial (cell array of matrix coefficents)
%
% OUTPUT
% E: eigenvalues of P
% V: eigenvectors of P associated with each eigval
%
% AUTHORS
% Written by Lucas Medina
% Minor modifications by Kayden Mimmack & Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets VPA digit precision higher for more accurate computations
digits(40);

% Generate symbolic linearization C of P
[C, ~] = Cp(P);
C = vpa(C);

% Get matrix V whose cols are right eigvecs, diagonal matrix D whose
% entries are corresponding eigvals 
[V, D] = eig(C);
V = vpa(V);
D = vpa(diag(D));

% Sort eigvals, eigvecs in order of increasing absolute value
[~, I] = sort(abs(D));
E = D(I);
V = V(:, I');
    