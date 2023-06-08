% var is an nx1 variable, where n is the number of grid cells
% M is the output map (nrows x ncols)
%
% it is possible the order of the lons/lats should be reversed in the
% output - see commented code

function [M, lons, lats] = make_map(basin, var)

[nrows, ncols] = size(basin.mask);

% bmask = (basin.mask);
bmask = flipud(basin.mask);

if sum(isnan(bmask(:))) == 0
    bmask = double(bmask);
    bmask(bmask==0) = NaN;
end

A = bmask(:);
inmask = ~isnan(A);

A(inmask) = var;
M = (reshape(A, nrows, ncols));
% M = flipud(reshape(A, nrows, ncols));

% figure,imagesc(M);

% tmp = (basin.grid_lat);
tmp = flipud(basin.grid_lat);
tmp = tmp(inmask);
lats = tmp;

% tmp = (basin.grid_lon);
tmp = flipud(basin.grid_lon);
tmp = tmp(inmask);
lons = tmp;

return
