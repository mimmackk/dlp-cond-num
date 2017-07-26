function Xm = XmComponent(m, P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% Constructs the Xm matrix described in the following paper:
% http://web.math.ucsb.edu/~mbueno/papers/GFPR.pdf
%
% INPUT
% m: the mth basis pencil to input matrix coeffs into
% P: cell array of matrix coefficients to plug into pencil
%
% OUTPUT
% Xm: A single Xm matrix used in constructing D_i basis
%
% AUTHORS
% Written by Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
k = length(P) - 1;
n = size(P{1}, 1);

%construction Um upper triangular matrix
u = 1 : k - m;
Um = sym(hankel(fliplr(u)));

%construction Lm lower triangular matrix
l = (k - m + 2) : (k + 1);
Lm = sym(fliplr(flipud(hankel(l))));

%instantiating Xm matrix frame
Xm = zeros(k, k);

%inserting lower triangular matrix in upper left quadrant of Xm
Xm((1 : m), (1 : m)) = Lm;

%adding rows/columns of 0's to Um to fit it into the lower right
%quadrant of Xm
for col = 1 : size(Xm((m + 1) : end, (m + 1) : end), 2) - size(Um, 2)
    %add columns to fix dimension
    Um = [Um zeros(size(Um, 1), 1)];
end

for col = 1 : size(Xm((m + 1) : end, (m + 1) : end), 1) - size(Um, 1)
    %add rows to fix dimension
    Um = [Um; zeros(1, size(Um, 2))];
end

%inserting Um into the lower right quadrant
Xm((m + 1) : end, (m + 1) : end) = Um;
Xm = sym(Xm);
replacementVector = sym(1 : (k + 1));

%replacing the frame reperesentations with block matrix
%coefficients from P
Xm = subs(Xm, replacementVector, P);

%multiplying the Um block matrices by -1
if k ~= m
    Xm((n * m + 1) : end, (n * m + 1) : end) = -1 * Xm((n * m + 1) : end, (n * m + 1) : end);    
end