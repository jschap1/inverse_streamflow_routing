% Calculate distance matrix
%
% Calculates the distance matrix for a set of points (x,y)

function d = calc_distance_matrix(x, y)

n = length(x);

% Calculate the distance matrix
d = zeros(n,n);
for i=1:n
    for j=1:n
%         d(i,j) = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);
        dist = distance(y(i), x(i), y(j), x(j)); % great circle distance
        d(i,j) = deg2km(dist);
    end
end

return