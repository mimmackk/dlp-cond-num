function main(epsilon, k, n, case_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
% main generates a random matrix polynomial of specified size and examines
% the condition number ratio of linearizations which are close to not being
% linearizations of P in one of four specified ways.
%
% INPUT
% epsilon: small value. The closer to 0 it is, the closer of a pencil in
%    DLP to a non-linearization we're examining.
% k: degree of our desired P
% n: size of our desired P
% case_num: value 1 - 4, specifies in what way we want D to be close to not
%    being a linearization for P
%
% AUTHORS
% Kayden Mimmack & Kaizen Towfiq, credits for graphing code to Lucas Medina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: add header docs to xmcomponent,
% upload all to github

% Sets digit precision higher for more accurate computations
digits(40);

[P, D] = generate_experiments(epsilon, k, n, case_num);

disp('Calculating condition number ratios for D in DLP, P...');
[ratios, eigvals] = cond_num_ratios(k - 1, D, P);
ratios = vpa(ratios);

% Graph our results
close all

% disp('RATIOS');
% disp(ratios);

% Index evenly spaced eigenvalues
eigvals_index = 1 : (n * k);
semilogy(eigvals_index, ratios, '-o');
hold on;

% legend('Condition Number Ratio D / P','location', 'best');
title('Condition Number Ratios \kappa_D(\delta) / \kappa_P(\delta)');

disp('RESULTS:');
disp(['k = ', num2str(k)]);
disp(['n = ', num2str(n)]);

% disp('Largest |eigval| = ');
% disp(vpa(abs(eigvals(n * k))));
% disp('Smallest |eigval| = ');
% disp(vpa(abs(eigvals(1))));
disp(double(eigvals));

% Calculates norm of each A_i
% for j = 0 : k
%     disp(['||A_', num2str(j), '|| = ']);
%     disp(norm(vpa(P{j + 1})));
% end

% Displays each A_i
for j = 0 : k
    disp(['||A_', num2str(j), '|| = ']);
    disp((double(P{j + 1})));
end


