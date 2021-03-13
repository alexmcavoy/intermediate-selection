function fp = fixation_probability(initial_state, final_state, states, transition_matrix)
% FIXATION_PROBABILITY
%   FIXATION_PROBABILITY(initial_state, final_state, transition_matrix)
%   takes as input an initial state, a final state, and a transition matrix
%   for an absorbing Markov chain and returns the probability of eventually
%   reaching the final (absorbing) state starting from the initial state.

    state_count = size(transition_matrix, 1);
    
    initial_index = get_index(initial_state, states);
    final_index = get_index(final_state, states);
    
    final_vector = zeros(state_count, 1);
    final_vector(final_index) = 1;
        
    B = transition_matrix-eye(state_count);
    B(1, 1) = 1;
    B(state_count, state_count) = 1;
    A = B;
    A(:, initial_index) = final_vector;
    
    fp = real(exp(logdet(A)-logdet(B)));
end