%% verify_cheb_matrix.m
%
% Verifies the Chebyshev differentiation matrix construction against known
% small cases (D_1, D_2) and tests the negative sum trick.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

fprintf('==============================================================\n');
fprintf('Chebyshev Differentiation Matrix Verification\n');
fprintf('==============================================================\n\n');

%% Verify small matrices

fprintf('1. Verifying small matrices (D_1, D_2)...\n');

% D_1: 2x2 matrix
[D1, x1] = cheb_matrix(1);
D1_exact = [0.5, -0.5; 0.5, -0.5];

if max(abs(D1(:) - D1_exact(:))) < 1e-14
    fprintf('   D_1 verification PASSED\n');
else
    fprintf('   D_1 verification FAILED\n');
    fprintf('   Computed:\n');
    disp(D1);
    fprintf('   Expected:\n');
    disp(D1_exact);
end

% D_2: 3x3 matrix
[D2, x2] = cheb_matrix(2);
D2_exact = [1.5, -2.0, 0.5; 0.5, 0.0, -0.5; -0.5, 2.0, -1.5];

if max(abs(D2(:) - D2_exact(:))) < 1e-14
    fprintf('   D_2 verification PASSED\n');
else
    fprintf('   D_2 verification FAILED\n');
    fprintf('   Computed:\n');
    disp(D2);
    fprintf('   Expected:\n');
    disp(D2_exact);
end

%% Verify negative sum trick

fprintf('\n2. Verifying negative sum trick (D * ones = zeros)...\n');

for N = [4, 8, 16, 32]
    [D, x] = cheb_matrix(N);
    result = D * ones(N+1, 1);
    max_error = max(abs(result));
    if max_error < 1e-12
        fprintf('   N = %d: PASSED (max error = %.2e)\n', N, max_error);
    else
        fprintf('   N = %d: FAILED (max error = %.2e)\n', N, max_error);
    end
end

%% Demonstrate differentiation accuracy

fprintf('\n3. Differentiation accuracy for u(x) = 1/(1+4x^2):\n');
fprintf('------------------------------------------------------\n');
fprintf('     N       Max Error\n');
fprintf('------------------------------------------------------\n');

u_func = @(x) 1.0 ./ (1.0 + 4.0 * x.^2);
u_prime = @(x) -8.0 * x ./ (1.0 + 4.0 * x.^2).^2;

for N = [4, 8, 16, 32, 64]
    [D, x] = cheb_matrix(N);
    v = u_func(x);
    w = D * v;
    w_exact = u_prime(x);
    error = max(abs(w - w_exact));
    fprintf('%6d %14.6e\n', N, error);
end

fprintf('------------------------------------------------------\n');

fprintf('\n==============================================================\n');
fprintf('Verification complete!\n');
fprintf('==============================================================\n');
