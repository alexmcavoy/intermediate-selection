function transition_matrix = build_matrix_pp(N, b, c, delta)
% BUILD_MATRIX_PP
%   BUILD_MATRIX_PP(N, b, c, delta) takes as input
%   the population size, N, the benefit of the good, b, the cost of the
%   good, c, and the selection intensity, delta. The output is the
%   transition matrix for pp-goods with these parameters on the cycle.
    
    % state i means that there are i-1 producers
    transition_matrix = zeros(N+1, N+1);
    
    % one producer in the population
    %
    % lose a producer
    transition_matrix(2, 1) = 1/N;
    % gain a producer
    payD = 0;
    payC = -2*c;
    transition_matrix(2, 3) = 2*(1/N)*(1/(1+exp(delta*(payD-payC))));
    
    % two producers in the population
    %
    % lose a producer
    payD = b;
    payC = b-2*c;
    transition_matrix(3, 2) = 2*(1/N)*(1/(1+exp(delta*(payC-payD))));
    % gain a producer
    payD = 0;
    payC = b-2*c;
    transition_matrix(3, 4) = 2*(1/N)*(1/(1+exp(delta*(payD-payC))));
    
    % N-2 producers in the population
    %
    % lose a producer
    payD = b;
    payC = 2*(b-c);
    transition_matrix(N-1, N-2) = 2*(1/N)*(1/(1+exp(delta*(payC-payD))));
    % gain a producer
    payD = b;
    payC = b-2*c;
    transition_matrix(N-1, N) = 2*(1/N)*(1/(1+exp(delta*(payD-payC))));
    
    % N-1 producers in the population
    %
    % lose a producer
    payD = 2*b;
    payC = 2*(b-c);
    transition_matrix(N, N-1) = 2*(1/N)*(1/(1+exp(delta*(payC-payD))));
    % gain a producer
    transition_matrix(N, N+1) = 1/N;

    for i=4:N-2
        % lose a producer
        payD = b;
        payC = 2*(b-c);
        transition_matrix(i, i-1) = 2*(1/N)*(1/(1+exp(delta*(payC-payD))));
        % gain a producer
        payD = 0;
        payC = b-2*c;
        transition_matrix(i, i+1) = 2*(1/N)*(1/(1+exp(delta*(payD-payC))));
    end
    
    for i=1:N+1
        transition_matrix(i, i) = 1-sum(transition_matrix(i, :));
    end

end