function transition_matrix = build_matrix_cf(N, b, c, delta)
% BUILD_MATRIX_CF
%   BUILD_MATRIX_CF(N, b, c, intensity) takes as input
%   the population size, N, the benefit of the good, b, the cost of the
%   good, c, and the selection intensity, delta. The output is the
%   transition matrix for cf-goods with these parameters on the cycle. 

    % state i means that there are i-1 producers
    transition_matrix = zeros(N+1, N+1);
    
    % one producer in the population
    %
    % lose a producer
    transition_matrix(2, 1) = 1/N;
    % gain a producer
    payD = [0, 0];
    payC = [-c, -c];
    for j=1:length(payD)
        transition_matrix(2, 3) = transition_matrix(2, 3) + ...
            2*(1/N)*(1/length(payD))*(1/(1+exp(delta*(payD(j)-payC(j)))));
    end
    
    % two producers in the population
    %
    % lose a producer
    payD = [b, 0];
    payC = [-c, b-c];
    for j=1:length(payD)
        transition_matrix(3, 2) = transition_matrix(3, 2) + ...
            2*(1/N)*(1/length(payD))*(1/(1+exp(delta*(payC(j)-payD(j)))));
    end
    % gain a producer
    payD = [0, 0];
    payC = [b-c, -c];
    for j=1:length(payD)
        transition_matrix(3, 4) = transition_matrix(3, 4) + ...
            2*(1/N)*(1/length(payD))*(1/(1+exp(delta*(payD(j)-payC(j)))));
    end
    
    % N-2 producers in the population
    %
    % lose a producer
    payD = [0, 0, b, b];
    payC = [2*b-c, b-c, b-c, -c];
    for j=1:length(payD)
        transition_matrix(N-1, N-2) = transition_matrix(N-1, N-2) + ...
            2*(1/N)*(1/length(payD))*(1/(1+exp(delta*(payC(j)-payD(j)))));
    end
    % gain a producer
    payD = [b, 0, b, 0];
    payC = [b-c, b-c, -c, -c];
    for j=1:length(payD)
        transition_matrix(N-1, N) = transition_matrix(N-1, N) + ...
            2*(1/N)*(1/length(payD))*(1/(1+exp(delta*(payD(j)-payC(j)))));
    end
    
    % N-1 producers in the population
    %
    % lose a producer
    payD = [0, 0, b, b, b, b, 2*b, 2*b];
    payC = [2*b-c, b-c, b-c, -c, 2*b-c, b-c, b-c, -c];
    for j=1:length(payD)
        transition_matrix(N, N-1) = transition_matrix(N, N-1) + ...
            2*(1/N)*(1/length(payD))*(1/(1+exp(delta*(payC(j)-payD(j)))));
    end
    % gain a producer
    transition_matrix(N, N+1) = 1/N;

    for i=4:N-2
        % lose a producer
        payD = [b, b, 0, 0];
        payC = [-c, b-c, b-c, 2*b-c];
        for j=1:length(payD)
            transition_matrix(i, i-1) = transition_matrix(i, i-1) + ...
                2*(1/N)*(1/length(payD))*(1/(1+exp(delta*(payC(j)-payD(j)))));
        end
        % gain a producer
        payD = [0, 0];
        payC = [-c, b-c];
        for j=1:length(payD)
            transition_matrix(i, i+1) = transition_matrix(i, i+1) + ...
                2*(1/N)*(1/length(payD))*(1/(1+exp(delta*(payD(j)-payC(j)))));
        end
    end
    
    for i=1:N+1
        transition_matrix(i, i) = 1-sum(transition_matrix(i, :));
    end

end