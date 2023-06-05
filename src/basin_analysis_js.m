% do all the basin analysis
%
% Very finicky about input flow direction file
% Must not have any cells that don't contribute at the outlet
% Sometimes it doesn't work then, either. Grr. 
% Can try adjusting search_max parameter for larger basins

function basin = basin_analysis_js(flow_dir, georef, gauge_list, gage, snapflag, lon, lat)

% new input: coordinates
basin = struct();
basin.lonv = lon;
basin.latv = lat;
      
r_earth = 6371.0072 * 1000;  % earth radius [m] (Ming's value)
% r_earth = 6378.1*1000; % radius of earth (m) (Jacob's value)

% basin grid georeferences
xllcorner    = georef(1);
yllcorner    = georef(2);
cellsize     = georef(3);
nodata_value = georef(4);

[nrows ncols] = size(flow_dir);
flow_dir(flow_dir==nodata_value) = NaN;

%% routing basin analysis

% basin mask
basin_mask = flow_dir*0+1;
ncells = nansum(nansum(basin_mask));

% lat and lon grids
grid_lat = repmat(flipud((yllcorner+((1:nrows)-0.5)*cellsize)'), 1, ncols).*basin_mask;
grid_lon = repmat(xllcorner+((1:ncols)-0.5)*cellsize, nrows, 1).*basin_mask;
% make sure these are right...

% grid area
grid_area = cellsize/180*pi*r_earth*cellsize/180*pi*r_earth*cos(grid_lat/180*pi); % m

% lat and lon of 1st order downstream, for stream plotting purpose
ds1_lat = grid_lat;
ds1_lat(flow_dir==1|flow_dir==2|flow_dir==8) = ds1_lat(flow_dir==1|flow_dir==2|flow_dir==8) + cellsize;
ds1_lat(flow_dir==4|flow_dir==5|flow_dir==6) = ds1_lat(flow_dir==4|flow_dir==5|flow_dir==6) - cellsize;

ds1_lon = grid_lon;
ds1_lon(flow_dir==2|flow_dir==3|flow_dir==4) = ds1_lon(flow_dir==2|flow_dir==3|flow_dir==4) + cellsize;
ds1_lon(flow_dir==6|flow_dir==7|flow_dir==8) = ds1_lon(flow_dir==6|flow_dir==7|flow_dir==8) - cellsize;

lats = [reshape(grid_lat(basin_mask==1), ncells, 1) reshape(ds1_lat(basin_mask==1), ncells, 1)];
lons = [reshape(grid_lon(basin_mask==1), ncells, 1) reshape(ds1_lon(basin_mask==1), ncells, 1)];

% downstream flow distance, great circle distance using haversine formula
hslat = (1-cos((grid_lat-ds1_lat)/180*pi))/2;
hslon = (1-cos((grid_lon-ds1_lon)/180*pi))/2;
flow_dist = 2*r_earth*asin(sqrt(hslat + cos(grid_lat/180*pi).*cos(ds1_lat/180*pi).*hslon));

% downstream grid search at all orders
search_max = 500; % tuning parameter
% search_max = 200; % tuning parameter
ds_i = zeros(nrows, ncols, search_max);
ds_j = zeros(nrows, ncols, search_max);
ds_length = zeros(nrows, ncols, search_max);

% initialize
ds_i(:, :, 1) = repmat((1:nrows)', 1, ncols);
ds_j(:, :, 1) = repmat(1:ncols, nrows, 1);
ds_flow_dir = flow_dir;
flow_dir_flat = flow_dir(:);
basin_mask_flat = basin_mask(:);

% flow length to outlet
flow_length = basin_mask*0;
flow_length_flat = flow_length(:);
flow_dist_flat = flow_dist(:);

% contributing area
grid_area_flat = grid_area(:);
flow_accum = grid_area;
flow_accum_flat = flow_accum(:);

for d=2:search_max
    
    ds_i_tmp = squeeze(ds_i(:, :, d-1));
    ds_j_tmp = squeeze(ds_j(:, :, d-1));
    
    ds_i_flat = ds_i_tmp(:);
    ds_j_flat = ds_j_tmp(:);
    
    % accumulate flow length
    flow_length_flat = flow_length_flat + flow_dist_flat(ds_i_flat+(ds_j_flat-1)*nrows);
    ds_length(:, :, d) = reshape(flow_length_flat, nrows, ncols);
    
    % move downstream
    ds_i_tmp(ds_flow_dir==1|ds_flow_dir==2|ds_flow_dir==8) = ds_i_tmp(ds_flow_dir==1|ds_flow_dir==2|ds_flow_dir==8) - 1;
    ds_i_tmp(ds_flow_dir==4|ds_flow_dir==5|ds_flow_dir==6) = ds_i_tmp(ds_flow_dir==4|ds_flow_dir==5|ds_flow_dir==6) + 1;
    
    ds_j_tmp(ds_flow_dir==2|ds_flow_dir==3|ds_flow_dir==4) = ds_j_tmp(ds_flow_dir==2|ds_flow_dir==3|ds_flow_dir==4) + 1;
    ds_j_tmp(ds_flow_dir==6|ds_flow_dir==7|ds_flow_dir==8) = ds_j_tmp(ds_flow_dir==6|ds_flow_dir==7|ds_flow_dir==8) - 1;
    
    ds_i(:, :, d) = ds_i_tmp;
    ds_j(:, :, d) = ds_j_tmp;
    
    ds_i_tmp_flat = ds_i_tmp(:);
    ds_j_tmp_flat = ds_j_tmp(:);
    
    % key step -- find the flow direction of the downstream cell
    try
        ds_flow_dir_flat = flow_dir_flat(ds_i_tmp_flat+(ds_j_tmp_flat-1)*nrows);
    catch
        error('annoying bug. grrr. idk y :(((((((')
    end
     ds_flow_dir = reshape(ds_flow_dir_flat, nrows, ncols);
    % accumulate contributing area
    % too bad the following doesn't work:
    % flow_accum_flat(ds_i_flat+(ds_j_flat-1)*nrows) = flow_accum_flat(ds_i_flat+(ds_j_flat-1)*nrows) + grid_area_flat;
    for i=1:numel(ds_i_flat)
        if ( ~isnan(ds_i_flat(i)) && (ds_i_flat(i)~=ds_i_tmp_flat(i) || ds_j_flat(i)~=ds_j_tmp_flat(i)) )
            %Y flow_accum_flat(ds_i_flat(i)+(ds_j_flat(i)-1)*nrows) = flow_accum_flat(ds_i_flat(i)+(ds_j_flat(i)-1)*nrows) + grid_area_flat(i);
            flow_accum_flat(ds_i_tmp_flat(i)+(ds_j_tmp_flat(i)-1)*nrows) = flow_accum_flat(ds_i_tmp_flat(i)+(ds_j_tmp_flat(i)-1)*nrows) + grid_area_flat(i);
        end
    end

    % all reach the outlet?
    if (nanstd(ds_i_tmp_flat.*basin_mask_flat)==0 && nanstd(ds_j_tmp_flat.*basin_mask_flat)==0)
        dmax = d;
        outlet_i = nanmean(ds_i_tmp_flat.*basin_mask_flat);
        outlet_j = nanmean(ds_j_tmp_flat.*basin_mask_flat);
        break;
    end
    
end

flow_length = reshape(flow_length_flat, nrows, ncols);
flow_accum = reshape(flow_accum_flat, nrows, ncols);

% delete extra downstream map
if exist('dmax')
    ds_i(:, :, dmax+1:search_max) = [];
    ds_j(:, :, dmax+1:search_max) = [];
    ds_length(:, :, dmax+1:search_max) = [];
    basin.dmax = dmax;
else
    warning('dmax does not exist. Basin may have non-contributing cells. Fix this.')
end

%% extract sub-basins for gauge stations

gauge_lat = gauge_list(:, 3);
gauge_lon = gauge_list(:, 4) - 360;

% lat/lons for basin
% lonlim = [min(grid_lon(:)), max(grid_lon(:))];
% latlim = [min(grid_lat(:)), max(grid_lat(:))];
    
plotfig = 1;
if plotfig
    % plot a figure showing gauge locations and flow directions
    figure
    plotraster(basin.lonv, basin.latv, flipud(flow_accum), '')
%     plotraster(lonlim, latlim, flipud(flow_accum), '')
    hold on
    plot(gauge_lon, gauge_lat, 'r.', 'MarkerSize', 15)
end

% Snap gauge locations to the flow direction-derived grid
if snapflag
    [gauge_lon, gauge_lat, gage] = snap_gauge_locations(gauge_lon, gauge_lat, grid_lon, grid_lat, gage, cellsize);
end

gauge_i = nrows-floor((gauge_lat-yllcorner)/cellsize);
gauge_j = floor((gauge_lon-xllcorner)/cellsize)+1;

if gauge_j == 0
    gauge_j = 1;
    % helps if the gauge location is off by one cell
    warning('Setting gauge_j to 1')
end

figure
plotraster(basin.lonv, basin.latv, flipud(flow_accum), 'facc')
hold on
plot(gauge_lon, gauge_lat, 'r.', 'markersize', 20)

% facc = flipud(flow_accum);
% ngauges = length(gauge_lon);
% for ii=1:ngauges
%     facc(gauge_i(ii), gauge_j(ii));
% end


gauge_accum = flow_accum(gauge_i+(gauge_j-1)*nrows);

% gauge_accum = flow_accum(logical(gauge_i+(gauge_j-1)*nrows));

%plot(gauge_area, gauge_accum, '.');

% remove stations with too inaccurate flow area
rm_gages = 0;
if rm_gages
    
    gauge_area = gauge_list(:, 5)*1609*1609; % convert from sq mi to sq m
    
    bad_gauges = gauge_area./gauge_accum>2 | gauge_area./gauge_accum<0.5;

    gauge_list(bad_gauges, :) = [];
    gauge_lat(bad_gauges) = [];
    gauge_lon(bad_gauges) = [];

    gauge_i(bad_gauges) = [];
    gauge_j(bad_gauges) = [];

    gauge_accum(bad_gauges) = [];
    gauge_area(bad_gauges) = [];
else
    gauge_area = gauge_accum; % convert from sq mi to sq m
end

%Y: remove duplicate gauges in one grid cell
gauge_i_j = [gauge_i,gauge_j];
[gauge_i_j_tmp, ia, ic] = unique(gauge_i_j,'rows','stable');
gauge_list = gauge_list(ia,:);
gauge_lat = gauge_lat(ia,:);
gauge_lon = gauge_lon(ia,:);
gauge_i = gauge_i(ia,:);
gauge_j = gauge_j(ia,:);
try
    gauge_accum = gauge_accum(ia,:);
    gauge_area = gauge_area(ia,:);
catch
    gauge_accum = gauge_accum(:,ia)';
    gauge_area = gauge_area(:,ia)';
end
   
% maybe we should keep downstream-most gauge? JRS

ngauges = numel(gauge_lat);

figure
plotraster(basin.lonv, basin.latv, flipud(flow_accum), 'Flow accumulation')
hold on
plot(gauge_lon, gauge_lat, 'r.', 'MarkerSize', 20)

% Snap to the flow accumulation grid ------------------------------

if snapflag
    
% Snap up to one cell away (thres = 1)
for ii=1:ngauges
    
    % check adjacent cells (VIC convention)
    f = zeros(9,1);
    try f(1) = flow_accum(gauge_i(ii)-1, gauge_j(ii)); end % T!
    try f(2)= flow_accum(gauge_i(ii)-1, gauge_j(ii)+1); end % TR!
    try f(3) = flow_accum(gauge_i(ii), gauge_j(ii)+1); end % R!
    try f(4) = flow_accum(gauge_i(ii)+1, gauge_j(ii)+1); end % BR
    try f(5) = flow_accum(gauge_i(ii)+1, gauge_j(ii)); end % B!
    try f(6) = flow_accum(gauge_i(ii)+1, gauge_j(ii)-1); end % BL
    try f(7) = flow_accum(gauge_i(ii), gauge_j(ii)-1); end % L!
    try f(8) = flow_accum(gauge_i(ii)-1, gauge_j(ii)-1); end % TL!
    f(9) = flow_accum(gauge_i(ii), gauge_j(ii));
    [fmax, maxind] = nanmax(f);
    
    switch maxind
        case 1 % T
            gauge_lon(ii) = gauge_lon(ii);
            gauge_lat(ii) = gauge_lat(ii) + cellsize;
        case 2 % TR
            gauge_lon(ii) = gauge_lon(ii) + cellsize;
            gauge_lat(ii) = gauge_lat(ii) + cellsize;            
        case 3 % R
            gauge_lon(ii) = gauge_lon(ii) + cellsize;
            gauge_lat(ii) = gauge_lat(ii);            
        case 4 % BR
            gauge_lon(ii) = gauge_lon(ii) + cellsize;
            gauge_lat(ii) = gauge_lat(ii) - cellsize;            
        case 5 % B
            gauge_lon(ii) = gauge_lon(ii);
            gauge_lat(ii) = gauge_lat(ii) - cellsize;            
        case 6 % BL
            gauge_lon(ii) = gauge_lon(ii) - cellsize;
            gauge_lat(ii) = gauge_lat(ii) - cellsize;            
        case 7 % L
            gauge_lon(ii) = gauge_lon(ii) - cellsize;
            gauge_lat(ii) = gauge_lat(ii);            
        case 8 % TL
            gauge_lon(ii) = gauge_lon(ii) - cellsize;
            gauge_lat(ii) = gauge_lat(ii) + cellsize;
    end
    
end

end

% remove duplicate gauges, keeping the gage with the largest flow
[gauge_lon, gauge_lat, gage, ia] = rm_duplicate_gauges(gauge_lon, gauge_lat, gage);
gauge_list = gauge_list(ia,:);
gauge_i = gauge_i(ia,:);
gauge_j = gauge_j(ia,:);

figure
plotraster(basin.lonv, basin.latv, (flipud(flow_accum)), 'Flow accumulation (snapped)')
hold on
plot(gauge_lon, gauge_lat, 'r.', 'MarkerSize', 20)

% gauge_accum = flow_accum(gauge_i+(gauge_j-1)*nrows);

gauge_accum = gauge_accum(ia,:);
gauge_area = gauge_area(ia,:);
ngauges = length(gauge_lon);

% sub-basin flow length

basin_mask_gauge = zeros(nrows, ncols, ngauges);
flow_length_gauge = zeros(nrows, ncols, ngauges);

for g=1:ngauges
    basin_mask_gauge(:, :, g) = basin_mask*0;
end

for i=1:nrows
    for j=1:ncols
        
        if (isnan(basin_mask(i, j)))
            continue;
        end
        
        ds_i_cell = squeeze(ds_i(i, j, :));
        ds_j_cell = squeeze(ds_j(i, j, :));
        
        for g=1:ngauges
            if (~isempty(find(ds_i_cell==gauge_i(g)&ds_j_cell==gauge_j(g),1)))
                basin_mask_gauge(i, j, g) = 1;
            end
        end
        
    end
end

for g=1:ngauges
    % if error, make sure gauges are snapped correctly
    flow_length_gauge(:, :, g) = (flow_length - flow_length(gauge_i(g), gauge_j(g))) .* squeeze(basin_mask_gauge(:, :, g));
end

%% Remove gauges that have NaN flow accumulation

% This happens due to a bug in the code, but it is only 3 out of 122
% gauges, so I am going to move on for now.

rm_ind = find(isnan(gauge_accum));

gauge_lon(rm_ind) = [];
gauge_lat(rm_ind) = [];
gauge_i(rm_ind) = [];
gauge_j(rm_ind) = [];
gauge_area(rm_ind) = [];
gauge_accum(rm_ind) = [];
basin_mask_gauge(:,:,rm_ind) = [];
gauge_list(rm_ind) = [];
flow_length_gauge(:,:,rm_ind) = [];
ia(rm_ind) = [];
gage(rm_ind,:) = [];

%% Write outputs to structure

basin.mask = basin_mask;
basin.grid_lat = grid_lat;
basin.grid_lon = grid_lon;
basin.flow_length = flow_length;
basin.grid_area = grid_area;
basin.flow_accum = flow_accum;
basin.ds_i = ds_i;
basin.ds_j = ds_j;
basin.ds_length = ds_length;
basin.gage_lat = gauge_lat;
basin.gage_lon = gauge_lon;
basin.gage_i = gauge_i;
basin.gage_j = gauge_j;
basin.gage_area = gauge_area;
basin.gage_accum = gauge_accum;
basin.mask_gage = basin_mask_gauge;
basin.gage_list = gauge_list;
basin.flow_length_gage = flow_length_gauge;
basin.ia = ia; % gauges to keep
basin.gage = gage; % discharge measurements

%% make lat/lon vectors

% This somehow gets the longitudes wrong (shifted)

% basin.lonv = reshape(basin.grid_lon, nrows*ncols, 1);
% basin.lonv(isnan(basin.lonv)) = [];

% basin.latv = reshape(basin.grid_lat, nrows*ncols, 1);
% basin.latv(isnan(basin.latv)) = [];   

      
return