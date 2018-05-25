function [classified, k, dist, idx] = fastKNN(trained, unknown, k, distance)

    if (nargin <= 3)
        % Minkowski Distance
        % for p == 2, Minkowski becomes equal Euclidean
        % for p == 1, Minkowski becomes equal city block metric
        % for p ~= 1 && p ~= 2 -> Minkowski https://en.wikipedia.org/wiki/Minkowski_distance
        distance = 2;
    end

    % trained data has one more column as unknown, the class
    [dist, idx] = getDistance(trained, unknown, distance);

    if (nargin <= 2)
        % determine k value when no one is given
        % possible number of categories + 1
        k = numel( unique( trained(:,end) ) ) + 1;
    end
    
    classified  = mode(trained(idx(1:k), end));


end % function fastKNN

function [values, idx] = getDistance(x, y, p)
    [values, idx] = sort (sum ( abs( (x(:,1:end-1) - y(ones(size(x,1),1), :) ) .^p), 2) .^ (1 / p));
end
