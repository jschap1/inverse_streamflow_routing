% Snaps gauge locations
%
% Snaps gauges to grid.
% Keeps gauge with largest flows in that grid cell

function [gauge_lon_new, gauge_lat_new, gage] = snap_gauge_locations(gauge_lon, gauge_lat, grid_lon, grid_lat, gage, cellsize)

grid_lon(isnan(grid_lon)) = 0;
lon = unique(grid_lon);
lon(end) = []; % hard-coded for 1/4 UMRB

grid_lat(isnan(grid_lat)) = 0;
lat = unique(grid_lat);
lat(1) = []; % hard-coded for 1/4 UMRB

res = cellsize;

m = length(gauge_lon); % number of gauges

gauge_lon_new = zeros(m,1);
gauge_lat_new = zeros(m,1);

% Snap the m gauge locations to the nearest grid cell center
% If a gauge is out of bounds of the study area, exclude it
for k=1:m
    
%     ismember([gauge_lon(k), gauge_lat(k)],grid_coords, 'rows')
    
%     find(min(abs(gauge_lon(k) - lon(:)))<(res/2));
%     find(min(abs(gauge_lat(k) - lat(:)))<(res/2));
    
    [val1, loc1] = min(abs(gauge_lon(k) - lon(:)));
    [val2, loc2] = min(abs(gauge_lat(k) - lat(:)));
    
    gauge_lon_new(k) = lon(loc1);
    gauge_lat_new(k) = lat(loc2);
    
%     gauge_lon_new(k) = lon(loc1) - res/2 + res;
%     gauge_lat_new(k) = lat(loc2) - res/2 - res;
    
    if val1 > res/2 || val2 > res/2
        warning(['gauge ' num2str(k) ' is out of bounds'])
        gauge_lon_new(k) = NaN;
        gauge_lat_new(k) = NaN;
    end
      
end

% Remove out of bounds gauges
rm_ind = unique([find(isnan(gauge_lon_new)), find(isnan(gauge_lat_new))]);

gauge_lon_new(rm_ind) = [];
gauge_lat_new(rm_ind) = [];
gage(rm_ind,:) = [];

% If there is more than one gauge in a grid cell, keep only the
% downstream-most gauge. This is calculated by assuming the gauge with the
% highest cumulative flow is the downstream-most gauge. 
%
% Note: this assumption could introduce some error, if the gauge data 
% are inaccurate.

[gauge_lon_new, gauge_lat_new, gage] = rm_duplicate_gauges(gauge_lon_new, gauge_lat_new, gage);

return