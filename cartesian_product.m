function r = cartesian_product(p, q)
% CARTESIAN_PRODUCT
%   CARTESIAN_PRODUCT(p, q) returns the product of a matrix, p, and a
%   column vector, q, where the rows are thought of as ordered tuples.

    [m, n] = size(p);
    r = zeros(m*length(q), n+1);
    row = 1;
    for i=1:m
        for j=1:length(q)
            r(row, 1:n) = p(i, :);
            r(row, n+1) = q(j);
            row = row+1;
        end
    end
end