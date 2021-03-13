function index = get_index(state, states)
% GET_INDEX
%   GET_INDEX(state, states) takes a state, represented by a row vector,
%   together with a matrix of states, in which each row represents a unique
%   state, and returns the row index of the state (and -1 if not found).

    index = -1;
    found = 0;
    current_index = 0;
    while current_index<size(states, 1) && ~found
        current_index = current_index + 1;
        if prod(state==states(current_index, :))==1
            index = current_index;
            found = 1;
        end
    end
end
