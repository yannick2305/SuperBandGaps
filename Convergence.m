%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [June 2025]
    Description:  [Gap Functions Convergence rates]
    --------------------------------------------------------------
%}

% --- Set the parameters ---

    Chain_length = 2000;
    N_values = fibonacci_substitution_lengths(Chain_length); 

    % --- Fix Lambda in the super-band gap
    lambda = 1;  

D_N_values = zeros(1, length(N_values));  % Will store D_N for each N

for idx = 1:length(N_values)
    N = N_values(idx);

    % Generate Fibonacci sequence
    sequence = generate_fibonacci_sequence(N);
    s = zeros(1, N);
    for i = 1:N
        s(i) = (sequence(i) == 'A') * l1 + (sequence(i) == 'B') * l2;
    end

    a = zeros(1, N);
    b = zeros(1, N);
    for n = 1:N
        s_prev = s(mod(n - 2, N) + 1);  % s_0 = s_N (periodic boundary)
        a(n) = 1/s_prev + 1/s(n);
        b(n) = -1/s(n);
    end

    D = zeros(1, N);
    D(1) = 0;
    D(2) = 1;

    for n = 3:N
        D(n) = (a(n-1) - lambda) * D(n-1) - b(n-2) * b(n-1) * D(n-2);
    end

    A = (-1)^N * prod(b);

    g_lambda = tridiag_det(a, b, lambda) - b(N)*b(N) * D(end);

    D_N_values(idx) = acosh( - g_lambda / (2*A) ) / (N + sum(s));
end

% Plot D_N_value as function of N
figure;
plot(N_values, D_N_values, 'r-o', 'LineWidth', 1.5);
xlabel('$N$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('$\beta_N(\lambda = 2)$', 'FontSize', 18, 'Interpreter', 'latex');
title('Lyapunov Exponent Approximation vs. $N$', 'FontSize', 20, 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
grid on;

% Compute successive increments (absolute differences)
increments = abs(diff(D_N_values));

% Plot the increments on a semilogy scale
figure;
semilogy(1:length(increments), increments, 'k-o', 'LineWidth', 1.5);
xlabel('Index $N$ of $\mathcal{F}_N$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('$|\beta_N(\lambda) - \beta_{N-1}(\lambda)|$', 'FontSize', 18, 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
grid on;


%%

function seq = generate_fibonacci_sequence(N)
    % Generate Fibonacci substitution sequence: A → AB, B → A
    seq = 'A';
    while length(seq) < N
        seq = regexprep(seq, 'A', 'X'); % temp marker
        seq = regexprep(seq, 'B', 'A');
        seq = regexprep(seq, 'X', 'AB');
    end
    seq = seq(1:N);
end


function dN = tridiag_det(a, b, lambda)

    % In the paper this is det(A_0-lambda)

    % Recursion formula from Wikipedia for tridiagonal matrix

    n = length(a);
    d = zeros(n+1, 1);
    d(1) = 1;             % i.e. d_0 = 1 
    d(2) = a(1) - lambda;
    
    for k = 3:n+1
        d(k) = (a(k-1) - lambda)*d(k-1) - b(k-2)^2 * d(k-2);
    end
    
    dN = d(n+1);
end


function lengths = fibonacci_substitution_lengths(N)
    % Initialize lengths array with the initial length (1 character: 'A')
    lengths = 1;
    
    % Initial counts: countA = 1 ('A'), countB = 0
    countA = 1;
    countB = 0;

    % Continue substitution until the total length reaches or exceeds N
    while lengths(end) < N
        % Apply substitution rules:
        % A → AB ⇒ countA_new = countA + countB (all A's become AB)
        % B → A  ⇒ countB_new = countA (all B's become A)
        newA = countA + countB;
        newB = countA;

        % Update counts
        countA = newA;
        countB = newB;

        % Record total length after substitution
        lengths(end+1) = countA + countB;
    end
end


