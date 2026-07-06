% MATLAB Programming Project - Spring 2025-2026
% Group Members:
% [Student 1 Danyah Adel - U23102095]
% [Student 2 Hadir Salah - U23100573]
% [Student 3 Hafsa Rahmanullah - U22100835]
% [Student 4 Salma Mohammed - U23200317]
%
% Description: This script performs two main tasks:
%   Task I: Solve a 5x5 linear system from an electrical circuit using:
%           - Gauss-Jordan Elimination (Exact Method)
%           - Jacobi Method (Iterative Method)
%           - Successive Over-Relaxation (SOR)
%           - Successive Under-Relaxation (SUR)
%           It also calculates the Condition Number (K) 
%           and Significant Digits (d).
%   Task II: Find root of a cubic polynomial using:
%            - False Position Method
%            - Secant Method

clear; clc; close all;

%% ========================================================================
%                                TASK I
%        Solve System of Linear Equations from an Electrical Circuit
%% ========================================================================

fprintf('--------------------------------------------------------------\n');
fprintf('            MATLAB Project - Spring 2025-2026                 \n');
fprintf('--------------------------------------------------------------\n\n');

% ========================================================================
% Define the System of Equations (Ax = B)
% From Kirchhoff''s Voltage Law (KVL) for the given circuit
% ========================================================================

% ========================================================================
% User Input for Resistances and Sources
% ========================================================================

fprintf('Enter the 13 resistance values in the correct order [R1 ... R13]:\n');
R = input('R = ');

fprintf('Enter the 5 source values in the correct order [V1 ... V5]:\n');
V = input('V = ');

% Check the correct number of values
if length(R) ~= 13 || length(V) ~= 5
    fprintf('\nInput size is incorrect. Resistances must be 13 and sources must be 5.\n');
    return;
end

% Check if all values are real
if ~isreal(R) || ~isreal(V)
    fprintf('\nAll input values must be real numbers.\n');
    return;
end

% Check if resistance values are positive
if any(R <= 0) || any(V <= 0)
    fprintf('All resistance and source values must be positive.\n');
    return;
end

% Build matrix A automatically from the circuit
A = [(R(1)+R(2)+R(3)+R(8)+R(9)), -R(3),       -R(9),       -R(8),        0;      % (Loop 1)
    -R(3), (R(3)+R(4)+R(5)+R(6)+R(7)),        0,          -R(7),       -R(6);    % (Loop 2)
    -R(9),        0,          (R(9)+R(10)+R(11)),         -R(11),       0;       % (Loop 3)
    -R(8),       -R(7),       -R(11), (R(8)+R(7)+R(11)+R(12)),         -R(12);   % (Loop 4)
     0,          -R(6),        0,          -R(12),      (R(6)+R(12)+R(13))];     % (Loop 5)

% Build vector B automatically from the source values
B = [V(1)-V(2);
     V(2)-V(3);
     V(4);
     0;
     -V(5)];

% Check if the system is singular before solving
% The determinant can be exactly zero or extremely close to zero
if det(A) == 0 || abs(det(A)) < 1e-12
    fprintf('The system is singular, so no unique solution exists.\n');
    return;
end

%% ------------------------------------------------------------------------
%  PART (A): Exact Solution using Gauss-Jordan Elimination
%% ------------------------------------------------------------------------
fprintf('\n--------------------------------------------------------------\n');
fprintf('TASK I - PART (A): Exact Solution using Gauss-Jordan Elimination\n');
fprintf('--------------------------------------------------------------\n');

% Gauss-Jordan Elimination
Aug = [A, B]; % Concatenation

% Get the size of the matrix
[n, ~] = size(A);

% Apply Gauss-Jordan elimination
for col = 1:n
     % Choose the best pivot row
    [~, pivot_row] = max(abs(Aug(col:n, col)));
    pivot_row = pivot_row + col - 1;

    % Swap rows if needed
    if pivot_row ~= col
        temp = Aug(col, :);
        Aug(col, :) = Aug(pivot_row, :);
        Aug(pivot_row, :) = temp;
    end

    % Make sure on-diagonal element is not zero
    if Aug(col, col) == 0
       fprintf('Zero on-diagonal element found, so the system cannot be solved using Gauss-Jordan Elimination Method.\n');
    return;
    end

    % Make the on-diagonal element 1
    % Divide the entire row by the current diagonal element
    Aug(col, :) = Aug(col, :) / Aug(col, col);

    % Make the rest of the column equal to zero
    for row = 1:n
        % Ignore the on-diagonal element
        if row ~= col
        % The value sitting above or below 1
            factor = Aug(row, col);
             % Make off-diagonal elements 0 for the specific row
            Aug(row, :) = Aug(row, :) - factor * Aug(col, :);
        end
    end
end

i_exact = Aug(:, end);

% Display the exact solution
fprintf('\nExact Solutions (i1 to i5):\n');
fprintf('i1 = %.10f A\n', i_exact(1));
fprintf('i2 = %.10f A\n', i_exact(2));
fprintf('i3 = %.10f A\n', i_exact(3));
fprintf('i4 = %.10f A\n', i_exact(4));
fprintf('i5 = %.10f A\n', i_exact(5));

%% ------------------------------------------------------------------------
%  PART (D): Condition Number Calculation
%% ------------------------------------------------------------------------
fprintf('\n--------------------------------------------------------------\n');
fprintf('TASK I - PART (D): Condition Number K\n');
fprintf('--------------------------------------------------------------\n');

% Calculate condition number using infinity norm
K = cond(A, inf);

fprintf('\nCondition Number K (inf-norm) = %.6e\n', K);

% System's health
threshold = input('Enter the threshold for which the system might be considered ill-conditioned: ');

if K <= threshold
    fprintf('System is WELL-CONDITIONED. Solutions are stable and reliable.\n');
else
    fprintf('System is ILL-CONDITIONED. Solutions are sensitive to small changes.\n');
end

%% ------------------------------------------------------------------------
%  PART (B): Iterative Methods
%% ------------------------------------------------------------------------
fprintf('\n--------------------------------------------------------------\n');
fprintf('TASK I - PART (B): Iterative Methods\n');
fprintf('--------------------------------------------------------------\n');

% Check Diagonal Dominance (Convergence Condition)
fprintf('\nChecking Convergence Condition (Diagonal Dominance):\n');

A_iter = A;
B_iter = B;

% At the beginning, assume the matrix is suitable for convergence
isDiagonallyDominant = true;
% Compare the diagonal element with the sum of the other elements
for i = 1:size(A_iter, 1)
    if abs(A_iter(i, i)) <= sum(abs(A_iter(i, :))) - abs(A_iter(i, i))
        isDiagonallyDominant = false;
        break;
    end
end

if isDiagonallyDominant
    fprintf('Matrix A is strictly diagonally dominant.\n');
    fprintf('Convergence is guaranteed for Jacobi, SOR, and SUR.\n');
else
    fprintf('Matrix A is NOT strictly diagonally dominant.\n');
    fprintf('Trying simple row reordering to improve convergence condition...\n');

    % Simple row reordering to improve diagonal entries
    for i = 1:size(A_iter,1)
        best_row = i;
        best_value = abs(A_iter(i,i));

        for r = i:size(A_iter,1)
            if abs(A_iter(r,i)) > best_value
                best_value = abs(A_iter(r,i));
                best_row = r;
            end
        end

        if best_row ~= i
            tempRow = A_iter(i,:);
            A_iter(i,:) = A_iter(best_row,:);
            A_iter(best_row,:) = tempRow;

            tempB = B_iter(i);
            B_iter(i) = B_iter(best_row);
            B_iter(best_row) = tempB;
        end
    end

    % Check again after reordering
    isDiagonallyDominant = true;
    for i = 1:size(A_iter, 1)
        if abs(A_iter(i, i)) <= sum(abs(A_iter(i, :))) - abs(A_iter(i, i))
            isDiagonallyDominant = false;
            break;
        end
    end

    if isDiagonallyDominant
        fprintf('\nAfter row reordering, the system became strictly diagonally dominant.\n');
    else
        fprintf('\nStrict diagonal dominance was still not fully achieved.\n');
        fprintf('We will continue with the reordered system and observe convergence numerically.\n');
    end
end

% Initial Conditions
x0 = [-0.05; -0.25; 0.5; -0.1; -0.9];  % i1, i2, i3, i4, i5

fprintf('\nInitial Conditions:\n');
fprintf('i1 = %.2f, i2 = %.2f, i3 = %.2f, i4 = %.2f, i5 = %.2f\n', ...
    x0(1), x0(2), x0(3), x0(4), x0(5));

% ========================================================================
% Jacobi Method
% ========================================================================

% Jacobi Method - 5 Iterations
fprintf('\n--------------------------------------------------------------\n');
fprintf('Jacobi Method (5 iterations):\n');
fprintf('--------------------------------------------------------------\n');

n = length(B_iter);
x = x0;
x_new = zeros(n, 1);

fprintf('\nIter\t  i1\t\t  i2\t\t  i3\t\t  i4\t\t  i5\t\t Step Size\n');
fprintf('___________________________________________________________________________________________________\n');

for k = 1:5
    x_old = x;
    for i = 1:n
        x_new(i) = (B_iter(i) - A_iter(i, [1:i-1, i+1:end]) * x_old([1:i-1, i+1:end])) / A_iter(i, i);
    end
    step_size = max(abs(x_new - x_old));
    fprintf('%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', k, x_new(1), x_new(2), x_new(3), x_new(4), x_new(5), step_size);
    x = x_new;
end
i_jacobi_5 = x_new;

% Jacobi Method - 20 Iterations
fprintf('\n--------------------------------------------------------------\n');
fprintf('Jacobi Method (20 iterations):\n');
fprintf('--------------------------------------------------------------\n');

x = x0;
x_new = zeros(n, 1);

fprintf('\nIter\t  i1\t\t  i2\t\t  i3\t\t  i4\t\t  i5\t\t Step Size\n');
fprintf('___________________________________________________________________________________________________\n');

for k = 1:20
    x_old = x;
    for i = 1:n
        x_new(i) = (B_iter(i) - A_iter(i, [1:i-1, i+1:end]) * x_old([1:i-1, i+1:end])) / A_iter(i, i);
    end
    step_size = max(abs(x_new - x_old));
    fprintf('%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', ...
        k, x_new(1), x_new(2), x_new(3), x_new(4), x_new(5), step_size);
    x = x_new;
end
i_jacobi_20 = x_new;

% ========================================================================
% Successive Over-Relaxation (SOR)
% ========================================================================
fprintf('\n--------------------------------------------------------------\n');
alpha_SOR = input('\nEnter alpha for SOR (1<alpha<2): ');
while alpha_SOR <= 1 || alpha_SOR >= 2
    fprintf('Invalid alpha. For SOR, alpha must be 1 < alpha < 2.\n');
    alpha_SOR = input('\nEnter alpha for SOR (1<alpha<2): \n');
end

fprintf('\n--------------------------------------------------------------\n');
fprintf('Successive Over-Relaxation (SOR) with alpha = %.2f\n', alpha_SOR);
fprintf('--------------------------------------------------------------\n');

% SOR - 5 iterations
fprintf('\nSOR (5 iterations):\n');

x = x0;

fprintf('\nIter\t  i1\t\t  i2\t\t  i3\t\t  i4\t\t  i5\t\t Step Size\n');
fprintf('___________________________________________________________________________________________________\n');

for k = 1:5
    x_old = x;
    for i = 1:n
        sum1 = A_iter(i, 1:i-1) * x(1:i-1);
        sum2 = A_iter(i, i+1:end) * x_old(i+1:end);
        x(i) = (1 - alpha_SOR) * x_old(i) + alpha_SOR * (B_iter(i) - sum1 - sum2) / A_iter(i, i);
    end
    step_size = max(abs(x - x_old));
    fprintf('%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', k, x(1), x(2), x(3), x(4), x(5), step_size);
end
i_sor_5 = x;

% SOR - 20 iterations
fprintf('\nSOR (20 iterations):\n');

x = x0;

fprintf('\nIter\t  i1\t\t  i2\t\t  i3\t\t  i4\t\t  i5\t\t Step Size\n');
fprintf('___________________________________________________________________________________________________\n');

for k = 1:20
    x_old = x;
    for i = 1:n
        sum1 = A_iter(i, 1:i-1) * x(1:i-1);
        sum2 = A_iter(i, i+1:end) * x_old(i+1:end);
        x(i) = (1 - alpha_SOR) * x_old(i) + alpha_SOR * (B_iter(i) - sum1 - sum2) / A_iter(i, i);
    end
    step_size = max(abs(x - x_old));
    fprintf('%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', k, x(1), x(2), x(3), x(4), x(5), step_size);
end
i_sor_20 = x;

% ========================================================================
% Successive Under-Relaxation (SUR)
% ========================================================================
fprintf('\n--------------------------------------------------------------\n');
alpha_SUR = input('\nEnter alpha for SUR (0<alpha<1): ');
while alpha_SUR <= 0 || alpha_SUR >= 1
    fprintf('Invalid alpha. For SUR, alpha must satisfy 0 < alpha < 1.\n');
    alpha_SUR = input('\nEnter alpha for SUR (0<alpha<1): \n');
end

fprintf('\n--------------------------------------------------------------\n');
fprintf('Successive Under-Relaxation (SUR) with alpha = %.2f\n', alpha_SUR);
fprintf('--------------------------------------------------------------\n');

% SUR - 5 iterations
fprintf('\nSUR (5 iterations):\n');

x = x0;

fprintf('\nIter\t  i1\t\t  i2\t\t  i3\t\t  i4\t\t  i5\t\t Step Size\n');
fprintf('___________________________________________________________________________________________________\n');

for k = 1:5
    x_old = x;
    for i = 1:n
        sum1 = A_iter(i, 1:i-1) * x(1:i-1);
        sum2 = A_iter(i, i+1:end) * x_old(i+1:end);
        x(i) = (1 - alpha_SUR) * x_old(i) + alpha_SUR * (B_iter(i) - sum1 - sum2) / A_iter(i, i);
    end
    step_size = max(abs(x - x_old));
    fprintf('%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', k, x(1), x(2), x(3), x(4), x(5), step_size);
end
i_sur_5 = x;

% SUR - 20 iterations
fprintf('\nSUR (20 iterations):\n');

x = x0;

fprintf('\nIter\t  i1\t\t  i2\t\t  i3\t\t  i4\t\t  i5\t\t Step Size\n');
fprintf('___________________________________________________________________________________________________\n');

for k = 1:20
    x_old = x;
    for i = 1:n
        sum1 = A_iter(i, 1:i-1) * x(1:i-1);
        sum2 = A_iter(i, i+1:end) * x_old(i+1:end);
        x(i) = (1 - alpha_SUR) * x_old(i) + alpha_SUR * (B_iter(i) - sum1 - sum2) / A_iter(i, i);
    end
    step_size = max(abs(x - x_old));
    fprintf('%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', k, x(1), x(2), x(3), x(4), x(5), step_size);
end
i_sur_20 = x;

%% ------------------------------------------------------------------------
%  PART (C): Significant Digits Calculation
%% ------------------------------------------------------------------------
fprintf('\n--------------------------------------------------------------\n');
fprintf('TASK I - PART (C): Significant Digits (d) for Approximations\n');
fprintf('--------------------------------------------------------------\n');

% Exact solution for the reordered system used in iterative methods
Aug_iter = [A_iter, B_iter];
[n_iter, ~] = size(A_iter);

for col = 1:n_iter
    [~, pivot_row] = max(abs(Aug_iter(col:n_iter, col)));
    pivot_row = pivot_row + col - 1;

    if pivot_row ~= col
        temp = Aug_iter(col, :);
        Aug_iter(col, :) = Aug_iter(pivot_row, :);
        Aug_iter(pivot_row, :) = temp;
    end

    if Aug_iter(col, col) == 0
        fprintf('Zero pivot encountered while solving the reordered system.\n');
        return;
    end

    Aug_iter(col, :) = Aug_iter(col, :) / Aug_iter(col, col);

    for row = 1:n_iter
        if row ~= col
            factor = Aug_iter(row, col);
            Aug_iter(row, :) = Aug_iter(row, :) - factor * Aug_iter(col, :);
        end
    end
end

i_exact_iter = Aug_iter(:, end);

% Jacobi - 5 iterations
fprintf('\nJacobi Method (5 iterations) - Significant Digits:\n');
for j = 1:5
    if i_exact_iter(j) == 0
        d = 0;
    else
        d_full = 1 - log10(2 * abs(i_exact_iter(j) - i_jacobi_5(j)) / abs(i_exact_iter(j)));
        if d_full == floor(d_full)
            d = d_full - 1;
        else
            d = floor(d_full);
        end
        if d < 0
            d = 0;
        end
    end
    fprintf('d%d = %.0f\n', j, d);
end

% Jacobi - 20 iterations
fprintf('\nJacobi Method (20 iterations) - Significant Digits:\n');
for j = 1:5
    if i_exact_iter(j) == 0
        d = 0;
    else
        d_full = 1 - log10(2 * abs(i_exact_iter(j) - i_jacobi_20(j)) / abs(i_exact_iter(j)));
        if d_full == floor(d_full)
            d = d_full - 1;
        else
            d = floor(d_full);
        end
        if d < 0
            d = 0;
        end
    end
    fprintf('d%d = %.0f\n', j, d);
end

% SOR - 5 iterations
fprintf('\nSOR (alpha = 1.2, 5 iterations) - Significant Digits:\n');
for j = 1:5
    if i_exact_iter(j) == 0
        d = 0;
    else
        d_full = 1 - log10(2 * abs(i_exact_iter(j) - i_sor_5(j)) / abs(i_exact_iter(j)));
        if d_full == floor(d_full)
            d = d_full - 1;
        else
            d = floor(d_full);
        end
        if d < 0
            d = 0;
        end
    end
    fprintf('d%d = %.0f\n', j, d);
end

% SOR - 20 iterations
fprintf('\nSOR (alpha = 1.2, 20 iterations) - Significant Digits:\n');
for j = 1:5
    if i_exact_iter(j) == 0
        d = 0;
    else
        d_full = 1 - log10(2 * abs(i_exact_iter(j) - i_sor_20(j)) / abs(i_exact_iter(j)));
        if d_full == floor(d_full)
            d = d_full - 1;
        else
            d = floor(d_full);
        end
        if d < 0
            d = 0;
        end
    end
    fprintf('d%d = %.0f\n', j, d);
end

% SUR - 5 iterations
fprintf('\nSUR (alpha = 0.8, 5 iterations) - Significant Digits:\n');
for j = 1:5
    if i_exact_iter(j) == 0
        d = 0;
    else
        d_full = 1 - log10(2 * abs(i_exact_iter(j) - i_sur_5(j)) / abs(i_exact_iter(j)));
        if d_full == floor(d_full)
            d = d_full - 1;
        else
            d = floor(d_full);
        end
        if d < 0
            d = 0;
        end
    end
    fprintf('d%d = %.0f\n', j, d);
end

% SUR - 20 iterations
fprintf('\nSUR (alpha = 0.8, 20 iterations) - Significant Digits:\n');
for j = 1:5
    if i_exact_iter(j) == 0
        d = 0;
    else
        d_full = 1 - log10(2 * abs(i_exact_iter(j) - i_sur_20(j)) / abs(i_exact_iter(j)));
        if d_full == floor(d_full)
            d = d_full - 1;
        else
            d = floor(d_full);
        end
        if d < 0
            d = 0;
        end
    end
    fprintf('d%d = %.0f\n', j, d);
end

fprintf('\n--------------------------------------------------------------\n');
fprintf('                         END OF TASK I                         \n');
fprintf('--------------------------------------------------------------\n');

%% ========================================================================
%                                TASK II
%              Root Finding Methods for Cubic Polynomial
%% ========================================================================

fprintf('--------------------------------------------------------------\n');
fprintf('                 TASK II: Root Finding Methods                \n');
fprintf('--------------------------------------------------------------\n\n');

maxIter = 100;

% Ask the user for the initial conditions
fprintf('Enter the interval [x0 x1]:\n');
initial_values = input('[x0 x1] = ');

% Check that exactly two values are entered
if length(initial_values) ~= 2
    fprintf('You must enter exactly two initial values.\n');
    return;
end

% Check that the values are real
if ~isreal(initial_values)
    fprintf('Initial values must be real numbers.\n');
    return;
end

x0 = initial_values(1);
x1 = initial_values(2);

% Tolerance
fprintf('Enter the tolerance:\n');
tol = input('Tolerance: ');

if ~isreal(tol)
    fprintf('Tolerance must be a real number.\n');
    return;
end

% Check if tolerance value is positive
if any(tol <= 0)
    fprintf('Tolerance must be positive.\n');
    return;
end

fprintf('\n');
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf(' TASK II: Root Finding Methods \n');
fprintf('══════════════════════════════════════════════════════════════════\n\n');

fprintf('Function: f(x) = (1/5)x^3 - (27/10)x^2 + 11x - 15\n');
fprintf('Interval: [%d, %d]\n', x0, x1);
fprintf('Tolerance: %.1f or less\n\n', tol);

%% ------------------------------------------------------------------------
% Part (A): False Position Method
%% ------------------------------------------------------------------------
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('False Position Method:\n');
fprintf('══════════════════════════════════════════════════════════════════\n');

iterations_fp = 0;
error_val = inf; % Initialized so the loop starts
brackets_fp = []; % Matrix that stores [lower, upper] limits

fprintf('\nIter\t x0\t\t x1\t\t xnew\t\t f(xnew)\t Step Size\n');
fprintf('__________________________________________________________________________________\n');

while (error_val > tol) && (iterations_fp < maxIter)
    iterations_fp = iterations_fp + 1;
    brackets_fp(iterations_fp, 1) = x0; 
    brackets_fp(iterations_fp, 2) = x1;

    f0 = equation(x0);
    f1 = equation(x1);

    if f0 * f1 >= 0
        fprintf('The interval does not help the program converge, enter another interval.\n')
        return;
    end

    xnew = x1 - f1 * ((x1 - x0) / (f1 - f0));
    fnew = equation(xnew);

    error_val = abs(fnew);

    % Display current iteration data
    fprintf('%2d\t %9.6f\t %9.6f\t %9.6f\t %9.6f\t %9.6f\n',iterations_fp, x0, x1, xnew, fnew, error_val);

    % Convergence Condition Check
    if f0 * fnew < 0
        x1 = xnew;
    else
        x0 = xnew;
    end
end

root_fp = xnew;
fprintf('\nTolerance reached after %d iterations.\n', iterations_fp);

fprintf('\nRoot found: x = %.6f\n', root_fp);
fprintf('Number of iterations: %d\n', iterations_fp);
fprintf('Final error: %.6f\n', abs(equation(root_fp)));

% Display stored brackets
fprintf('\nStored Brackets [Lower Limit, Upper Limit]:\n');
fprintf('Iter\t Lower\t\t Upper\n');
for i = 1:size(brackets_fp, 1)
fprintf('%2d\t %9.6f\t %9.6f\n', i, brackets_fp(i, 1), brackets_fp(i, 2));
end

%% ------------------------------------------------------------------------
% Part (B): Secant Method
%% ------------------------------------------------------------------------
fprintf('\n══════════════════════════════════════════════════════════════════\n');
fprintf('Secant Method:\n');
fprintf('══════════════════════════════════════════════════════════════════\n');

x0 = initial_values(1);
x1 = initial_values(2);

% Variables for tracking
iterations_sec = 0;
error_val = inf;      % Initialization so the loop starts
brackets_sec = [];    % Matrix that stores [x0, x1] pairs

fprintf('\nIter\t x0\t\t x1\t\t xnew\t\t f(xnew)\t Step Size\n');
fprintf('__________________________________________________________________________________\n');

while (error_val > tol) && (iterations_sec < maxIter)
    iterations_sec = iterations_sec + 1;
    
    brackets_sec(iterations_sec, 1) = x0; 
    brackets_sec(iterations_sec, 2) = x1;
    
    f0 = equation(x0);
    f1 = equation(x1);
    
    if f1 - f0 == 0
        fprintf('The given interval causes the denominator to become zero, enter another interval.\n')
        return;
    end

    xnew = x1 - ((x1 - x0) / (f1 - f0)) * f1;
    fnew = equation(xnew);
    
    error_val = abs(fnew);
    
    % Display current iteration data
    fprintf('%2d\t %9.6f\t %9.6f\t %9.6f\t %9.6f\t %9.6f\n', iterations_sec, x0, x1, xnew, fnew, error_val);
    
    x0 = x1;
    x1 = xnew;
end

root_sec = xnew;

fprintf('\nRoot found: x = %.6f\n', root_sec);
fprintf('Number of iterations: %d\n', iterations_sec);
fprintf('Final error: %.6f\n', abs(equation(root_sec)));

%% ------------------------------------------------------------------------
% Compare which of the 2 methods is quicker to find the root
%% ------------------------------------------------------------------------
fprintf('\n══════════════════════════════════════════════════════════════════\n');
fprintf('Comparison:\n');
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('False Position Method: %2d iterations\n', iterations_fp);
fprintf('Secant Method: %2d iterations\n', iterations_sec);
if iterations_sec < iterations_fp
fprintf('\nSecant Method is FASTER (fewer iterations).\n');
elseif iterations_fp < iterations_sec
fprintf('\nFalse Position Method is FASTER (fewer iterations).\n');
else
fprintf('\nBoth methods required the same number of iterations.\n');
end

function f = equation(x)
    f = (1/5)*x.^3 - (27/10)*x.^2 + 11*x - 15;
end
fprintf('\n--------------------------------------------------------------\n');
fprintf('                         END OF TASK II                         \n');
fprintf('--------------------------------------------------------------\n');