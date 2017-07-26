function D = DLPbasis(m, P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% DLPbasis returns the matrix coefficients for the mth DLP basis
% linearization.
%
% INPUT
% m: the mth basis pencil to input matrix coeffs into
% P: cell array of matrix coefficients to plug into pencil
%
% OUTPUT
% D: cell array consisting of the matrix coefficients for mth basis pencil
%    in DLP
%
% AUTHORS
% Written by Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
k = length(P) - 1;

if 1 <= m && m <= k
    X1 = XmComponent(m, P);
    X0 = XmComponent(m - 1, P);
    D = {X0, X1};
else
    fprintf('m-value is unacceptable, 1<=m<=k')
end