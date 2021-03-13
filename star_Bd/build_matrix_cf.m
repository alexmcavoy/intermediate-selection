function transition_matrix = build_matrix_cf(N, b, c, intensity)
% BUILD_MATRIX_CF
%   BUILD_MATRIX_CF(N, b, c, intensity) takes as input
%   the population size, N, the benefit of the good, b, the cost of the
%   good, c, and the selection intensity. The output is the transition
%   matrix for cf-goods with these parameters on the star.
    
    states = cartesian_product(transpose(0:1:1), transpose(0:1:N-1));
    state_count = size(states, 1);

    transition_matrix = zeros(state_count, state_count);
    for i=1:state_count
        current_state = states(i, :);
        for leaf_recipient=1:N-1
            if current_state(1)==1
                payoff_hub = current_state(2)*b - c;
                payoffC_leaf = -c;
                payoffD_leaf = 0;
                if leaf_recipient<=current_state(2)
                    payoff_leaf_recipient = b-c;
                    F = exp(intensity*([payoff_hub, payoffC_leaf, payoffD_leaf, payoff_leaf_recipient]-...
                        max([payoff_hub, payoffC_leaf, payoffD_leaf, payoff_leaf_recipient])));
                    F(2) = (current_state(2)-1)*F(2);
                    F(3) = (N-1-current_state(2))*F(3);
                else
                    payoff_leaf_recipient = b;
                    F = exp(intensity*([payoff_hub, payoffC_leaf, payoffD_leaf, payoff_leaf_recipient]-...
                        max([payoff_hub, payoffC_leaf, payoffD_leaf, payoff_leaf_recipient])));
                    F(2) = current_state(2)*F(2);
                    F(3) = (N-2-current_state(2))*F(3);
                end
                for j=1:state_count
                    next_state = states(j, :);
                    if sum(abs(current_state-next_state))==1
                        index = find(current_state~=next_state);
                        if index==1
                            if leaf_recipient<=current_state(2)
                                transition_matrix(i, j) = transition_matrix(i, j) + ...
                                    (1/(N-1))*(F(3)/sum(F));
                            else
                                transition_matrix(i, j) = transition_matrix(i, j) + ...
                                    (1/(N-1))*((F(3)+F(4))/sum(F));
                            end
                        else
                            if 1+current_state(2)==next_state(2)
                                transition_matrix(i, j) = transition_matrix(i, j) + ...
                                    (1/(N-1))*(F(1)/sum(F))*(1-current_state(2)/(N-1));
                            end
                        end
                    end
                end
            else
                payoff_hub = current_state(2)*b;
                payoffC_leaf = -c;
                payoffD_leaf = 0;
                F = exp(intensity*([payoff_hub, payoffC_leaf, payoffD_leaf]-...
                    max([payoff_hub, payoffC_leaf, payoffD_leaf])));
                F(2) = current_state(2)*F(2);
                F(3) = (N-1-current_state(2))*F(3);
                for j=1:state_count
                    next_state = states(j, :);
                    if sum(abs(current_state-next_state))==1
                        index = find(current_state~=next_state);
                        if index==1
                            transition_matrix(i, j) = transition_matrix(i, j) + ...
                                (1/(N-1))*(F(2)/sum(F));
                        else
                            if current_state(2)==1+next_state(2)
                                transition_matrix(i, j) = transition_matrix(i, j) + ...
                                    (1/(N-1))*(F(1)/sum(F))*(current_state(2)/(N-1));
                            end
                        end
                    end
                end
            end
        end
    end
    
    for i=1:state_count
        transition_matrix(i, i) = 1-sum(transition_matrix(i, :));
    end

end