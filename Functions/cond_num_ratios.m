function [ratios, eigvals] = cond_num_ratios(h, L, P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% Finds ratio of condition number of linearization L over cond num of 
% matrix polynomial P for every shared eigenvalue.
%
% INPUT
% h: indexing variable that indicates where in eigvec for L to extract
%    eigvec for P
% L: Linearization of P
% P: Matrix polynomial (cell array of matrix coefficents A_0 to A_k)
% family: String which specifies which type of linearization L is: DLP, Tp,
%    T, or Ds
%
% OUTPUT
% ratios: Array of condition number ratios kL / kP
% eigvals: Array of eigenvalues shared by L, P
%
% AUTHORS
% Written by Lucas Medina
% Modifications by Kayden Mimmack & Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets VPA digit precision higher for more accurate computations
digits(40);

% Convert L to VPA format, get degree and size of P
L = {vpa(L{1}), vpa(L{2})};
k = length(P) - 1;
n = size(P{1}, 1);

% Initializes array of ratios for each condition number
ratios = zeros(1, n * k);
ratios = sym(ratios, 'f');

% Linearize L into one single matrix C and get its eigenpairs
[C, L1_inverse] = Cp(L);
C = vpa(C);
[reigvecsL, eigvals] = eig(C);
reigvecsL = vpa(reigvecsL);

% Order eigenvalues & right eigenvectors of L by ascending modulus
eigvals = vpa(diag(eigvals));
[~, I] = sort(abs(eigvals));
eigvals = eigvals(I);
reigvecsL = reigvecsL(:, I');

% Get left eigenvectors of L
[leigvecsL, left_eigvals] = eig(C.');
left_eigvals = vpa(diag(left_eigvals));
I = zeros(1, length(eigvals));

% Find the equivalent ordering of left and right eigvecs for same eigvals
for j = 1 : length(eigvals)
    index = find(eval(left_eigvals) == eval(eigvals(j)), 1);
    I(j) = index;
end

% Rearrange left eigenvectors to match ordering of right eigenvectors
leigvecsL = vpa(leigvecsL);
leigvecsL = leigvecsL(:, I');
leigvecsL = ctranspose(L1_inverse) * conj(leigvecsL);

% Initialize matrices of eigenvectors of P
reigvecsP = zeros(n, n * k);
reigvecsP = sym(reigvecsP, 'f');
leigvecsP = reigvecsP;

% Extract eigenvectors of P from eigenvectors of L
for j = 1 : length(eigvals)
    
    start_index = n * (k - h - 1) + 1;
    end_index   = n * (k - h);
    
    reigvecsP(:, j) = reigvecsL(start_index : end_index, j);
    leigvecsP(:, j) = leigvecsL(start_index : end_index, j);

end

% Calculate condition number for each eigenvalue
for j = 1 : length(eigvals)

    curr_eigval = eigvals(j);
    
    % Get corresponding right, left eigvecs for L
    xL = reigvecsL(:, j);
    yL = leigvecsL(:, j);
    
    % Get corresponding right, left eigvecs for P
    xP = reigvecsP(:, j);
    yP = leigvecsP(:, j);

    % Calculate condition numbers for L, P
    kL = cond_num(L, curr_eigval, xL, yL);
    kP = cond_num(P, curr_eigval, xP, yP);

    disp(['Current eigenvalue: ', num2str(double(curr_eigval))]);
    disp(['Condition number of P: ', num2str(double(kP))]);
    disp(['Condition number of L: ', num2str(double(kL))]);

    ratios(j) = vpa(kL / kP);

    status = double(j / (n * k) * 100);
    disp([num2str(status), '% complete']);

end