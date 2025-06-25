%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [June 2025]
    Description:  [Gap Functions Convergence]
    --------------------------------------------------------------
%}

close all;
clear all;

% --- Parameters ---

Chain_length = 2000;

N_values = 1:16; %fibonacci_substitution_lengths(Chain_length); 
N_lam = 1000;
l1 = 1;
l2 = 2;
lambdas = linspace(0, 2.6, N_lam);

% --- Setup Video Writer (MP4) ---
v = VideoWriter('FibonacciTilingGapFunctionsFixed', 'MPEG-4');
v.FrameRate = 3;  % Adjust as needed
open(v);

figure('Position', [100 100 800 600]);

% --- Loop Over N Values ---
for idx = 1:length(N_values)
    N = N_values(idx);
    
    % Generate Fibonacci sequence
    sequence = generate_fibonacci_sequence(N);
    s = zeros(1, length(sequence));
    for i = 1:length(s)
        s(i) = (sequence(i) == 'A') * l1 + (sequence(i) == 'B') * l2;
    end
    disp(length(s));

    a = zeros(1, length(sequence));
    b = zeros(1, length(sequence));
    for n = 1:length(sequence)
        s_prev = s(mod(n - 2, length(sequence)) + 1);  % s_0 = s_N (periodic boundary)
        a(n) = 1/s_prev + 1/s(n);
        b(n) = -1/s(n);
    end

    D_N_values = zeros(size(lambdas));


    for k = 1:length(lambdas)
        lambda = lambdas(k);
        D = zeros(1, length(sequence));

        D(1) = 0;
        D(2) = 1;

        for n = 3:length(sequence)
            D(n) = (a(n-1) - lambda) * D(n-1) - b(n-2) * b(n-1) * D(n-2);
        end
    
        A = (-1)^length(sequence) * prod(b);

        g_lambda = tridiag_det(a, b, lambda) - b(length(sequence))*b(length(sequence)) * D(end);
        
        D_N_values(k) = acosh( - g_lambda / (2*A) ) / ( length(sequence) + sum(s) );
    end

    % Plot and capture frame
    clf;
    plot(lambdas, D_N_values, 'b', 'LineWidth', 1.5);
    xlabel('$\lambda$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$\beta_N(\lambda)$', 'FontSize', 18, 'Interpreter', 'latex');
    title(['Fibonacci recursion level: $N = ' num2str(N) '$'], 'FontSize', 20, 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    ylim([0, 0.15]);
    grid on;

    frame = getframe(gcf);
    writeVideo(v, frame);
end

% --- Finalize ---
close(v);
disp('Movie created: FibonacciGapDecay.avi');


%% --- Defining Functions ---


function seq = generate_fibonacci_sequence(N)
    % Generate Fibonacci substitution sequence with N substitution steps
    % Rule: A → AB, B → A
    seq = 'A';
    for k = 1:N
        % Apply substitution to each character in the current sequence
        next_seq = '';
        for i = 1:length(seq)
            if seq(i) == 'A'
                next_seq = [next_seq 'AB'];
            else  % seq(i) == 'B'
                next_seq = [next_seq 'A'];
            end
        end
        seq = next_seq;
    end
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

