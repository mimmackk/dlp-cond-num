function [P, D] = generate_experiments(epsilon, k, n, case_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% generate_experiments creates a random matrix polynomial P and a
% linearization D of P in DLP which satisfy the conditions for each of our
% cases that qualify D as being close to not a linearization
%
% INPUT
% epsilon: small value. The closer to 0 it is, the closer D is to not being
%   a linearization.
% k: degree of our desired P
% n: size of our desired P
% case_num: value 1 - 4, specifies in what way we want D to be close to not
%   being a linearization for P
%
% OUTPUT
% P: Matrix polynomial (cell array of matrix coefficents)
% D: Linearization of P in DLP epsilon-close to not being a linearization
%
% AUTHORS
% Kayden Mimmack & Kaizen Towfiq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets VPA digit precision higher for more accurate computations
digits(40);

% Case 1: Fixed P, construct v close to nearby v' such that D(P, v') is 
% not a linearization.
if case_num == 1
    
    P = rand_matrix_poly(k, n, 50);
    [eigvals, ~] = eigpairs(P);

    % Perterb median eigenvalue
    delta = eigvals(ceil(length(eigvals) / 2)) + epsilon;
    
    % Generate D_{k - 1} - delta * D_k, which is close to nonlinearization
    % bc delta is close to an eigenvalue of P
    D  = DLPbasis(k - 1, P);
    Dk = DLPbasis(k, P);
    D{1} = D{1} - delta * Dk{1};
    D{2} = D{2} - delta * Dk{2};
    
% Case 2: Fixed P with infinite eigenvalue, construct v such that v_1 is
% close to zero.
elseif case_num == 2
    

% Case 3: Fixed v, construct P close to nearby P' such that D(P', v) is
% not a linearization. Here using D = D_1, so want A_0 close to singular.
elseif case_num == 3
    
    P = rand_matrix_poly(k, n, 50); 
        
    % Fill last row, column of P with zeros so P is block diag
    for j = 1 : (k + 1)
        P{j}(:, n) = 0;
        P{j}(n, :) = 0;
    end

    % Put x^k - delta in lower right corner -> epsilon^(1 / k) is eigval
    P{k + 1}(n, n) = 1;
    P{1}(n, n) = -1 * epsilon;

    % Generate D1
    D = DLPbasis(1, P);
    D{1} = -1 * D{1};
    
% Case 4: Fixed v with v_1 = 0, construct P with very large eigenvalue.
% Here using D = D_k, so want A_k close to singular.
elseif case_num == 4
    
    P = rand_matrix_poly(k, n, 50); 
    
    % Copy 2nd to last row of A_k to last row, add epsilon to corner entry
    % Small epsilon -> A_k close to singular
    P{k + 1}(n, :) = P{k + 1}(n - 1, :);
    P{k + 1}(n, 1) = P{k + 1}(n, 1) + epsilon;

    % Generate Dk
    D = DLPbasis(k, P);
    D{1} = -1 * D{1};
    
end
    
    