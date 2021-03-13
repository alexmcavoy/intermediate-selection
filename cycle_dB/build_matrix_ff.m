function transition_matrix = build_matrix_ff(N, b, c, intensity)
% BUILD_MATRIX_FF
%   BUILD_MATRIX_FF(N, b, c, intensity) takes as input
%   the population size, N, the benefit of the good, b, the cost of the
%   good, c, and the selection intensity. The output is the transition
%   matrix for ff-goods with these parameters on the cycle
    
    transition_matrix = build_matrix_pp(N, b/2, c/2, intensity);

end