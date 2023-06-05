% Generate SWOT-like measurements for the Ohio basin
%
% 1/11/2023
% Likely will need to do parts of this in R or GRASS
%   -Generating random errors with CoSMoS
%   -Intersecting the SWOT orbit with the Ohio river network
%
% INPUTS
% River pixel coordinates
% SWOT orbit data
% GRWL river width data
% Ohio flow directions map OR basin mask
%
% Requires GDAL

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

fs = 18;
lw = 2;
ms = 20;

%% Load Ohio basin data

load('./ohio_data/setup-1-gage.mat')

%% Get river pixels

A = load('./ohio_data/river_pixel_coordinates.txt');
rlon = A(:,1); rlat = A(:,2);

clearvars Data
Data.Geometry  = 'Point';
Data.Name = 'centerline';
Data.Lon = rlon;
Data.Lat = rlat;

figure
geoshow(Data);

save('./ohio_data/river_pixels_1_16.mat', 'Data')

%% Get SWOT orbit

swotfile = './ohio_data/swot_ohio.shp';
swot_orbit = shaperead(swotfile);

% Crop to domain
% ogr2ogr -clipsrc -90 34 -76 43 ./basin_setup/swot_ohio.shp /Volumes/HD3/ISR/06_Data/Swaths/30s/sph_science_nadir/swot_science_orbit_sept2015-v2_nadir.shp

figure
mapshow(swot_orbit)

%% Get river widths

grwlfile = './ohio_data/GRWL_summaryStats_V01.01/GRWL_summaryStats.shp';
grwl = shaperead(grwlfile);

figure(2)
mapshow(grwl)
hold on
mapshow(Data)
title('Ohio River Pixels and GRWL')
xlabel('Lon')
ylabel('Lat')
set(gca, 'fontsize', fs)

% Crop GRWL to the extent of the Ohio river pixels
% Make a bounding box shapefile

bbfile = './ohio_data/boundingbox_for_cropping.shp';

x = [-90 -76 -76 -90]  ;
y = [34 34 42 42] ;
bb.Geometry = 'Polygon' ;
bb.X = x  ;  % latitude
bb.Y = y ;  % longitude
bb.Name = 'Rectangle' ; 

figure
mapshow(bb)

shapewrite(bb, bbfile)

% Crop in GDAL
!cd /Volumes/HD3/ISR/inverse_streamflow_routing
!ogr2ogr -clipsrc -90 34 -76 43 ./ohio_data/grwl_ohio.shp ./ohio_data/GRWL_summaryStats_V01.01/GRWL_summaryStats.shp
cropped_grwl = shaperead('./ohio_data/grwl_ohio.shp');
figure
mapshow(cropped_grwl)

% Get index for GRWL reaches > 100 m
% Get index for GRWL reaches > 50 m
nr = length(cropped_grwl); % number of reaches
id100 = zeros(nr, 1);
id50 = zeros(nr, 1);
for i=1:nr
    if cropped_grwl(i).width_mean > 100
        id100(i) = 1;
        id50(i) = 1;
    elseif cropped_grwl(i).width_mean > 50
        id50(i) = 1;
    end
end

figure
subplot(1,2,1)
mapshow(cropped_grwl(logical(id50)))
hold on
mapshow(Data)
title('Reaches wider than 50m')
subplot(1,2,2)

symbspec = makesymbolspec('Line', {'Default', 'Color', 'red', 'LineWidth', lw});

figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'Ohio')
hold on
mapshow(cropped_grwl(logical(id100)), 'SymbolSpec', symbspec)
mapshow(swot_orbit)
title('Reaches wider than 100m')

%% Find locations where SWOT crosses river centerlines

grwl50 = cropped_grwl(logical(id50));
grwl100 = cropped_grwl(logical(id100));

P = [[grwl100(:).X]; [grwl100(:).Y]]'; % "points to search" 
res = 0.25; % resolution (degrees)

intersection_points_all = zeros(0,3); % third col is day of overpass
for i=1:length(swot_orbit)

    query_points = [swot_orbit(i).X; swot_orbit(i).Y]';
    
    % need to densify the query points
    [latout,lonout] = interpm(query_points(:,1),query_points(:,2),res);
    query_points = [latout,lonout];
    
    [k,dist] = dsearchn(P,query_points);
    intersection_points = P(k,:);
    intersection_points = intersection_points(dist<res/2,:);
    
    intersection_day = str2double(swot_orbit(i).START_TIME(5:6));
    intday = repmat(intersection_day, size(intersection_points,1),1);
    
    intersection_points_all = [intersection_points_all; [intersection_points,intday]];
    
end

figure
mapshow(swot_orbit)
hold on
mapshow(grwl100, 'SymbolSpec', symbspec)
plot(intersection_points_all(:,1), intersection_points_all(:,2), 'blue.', 'markersize', ms)

dlmwrite('./ohio_data/swot_overpass_points_100m.txt',intersection_points_all)

T = table(intersection_points_all(:,1),intersection_points_all(:,2),intersection_points_all(:,3));
T.Properties.VariableNames{1} = 'lon';
T.Properties.VariableNames{2} = 'lat';
T.Properties.VariableNames{3} = 'orbit_day';

% Check if overpass location is within basin extent

fdfile = './ohio_data/ohio_flowdir.tif';
[fdir, R, lon, lat] = geotiffread2(fdfile);

ind = zeros(size(T,1),1);
j1 =  zeros(size(T,1),1);
j2 =  zeros(size(T,1),1);

for i=1:size(T,1)
    
    [j1(i),j2(i)] = latlon2pix(R, T{i,2}, T{i,1});
    j1(i) = round(j1(i));
    j2(i) = round(j2(i));
    
    try
    if isnan(basin.mask(j1(i),j2(i)))
        ind(i) = 0;
    else
        ind(i) = 1;
    end
    catch
        ind(i) = 0;
    end
    
end

ind = logical(ind);

figure
plotraster(basin.lonv, basin.latv, flipud(basin.mask), 'Virtual gauges (100 m)')
hold on
mapshow(cropped_grwl(logical(id100)), 'SymbolSpec', symbspec, 'color', 'blue')
mapshow(swot_orbit)
plot(intersection_points_all(ind,1),intersection_points_all(ind,2), '.r', 'markersize', ms)

T1 = [T,table(ind, j1, j2)];
T1.Properties.VariableNames{4} = 'in_basin';
T1.Properties.VariableNames{5} = 'row_ind';
T1.Properties.VariableNames{6} = 'col_ind';

T2 = T1(T1.in_basin==1,:);
m = size(T2,1); % number of virtual gauges

%% Create gauge list and route true runoff to each SWOT overpass location

nt = 365;
gauge_list = [ones(m,1), 1000+(1:m)', T2{:,2}, T2{:,1} + 360, 1000*ones(m,3)];

fdfile = './ohio_data/ohio_flowdir.tif';
[fdir, R, lon, lat] = geotiffread2(fdfile);
fdir = flipud(fdir);
figure,plotraster(lon, lat, fdir, 'flow directions');

timestep = 86400; % s
[nrows, ncols] = size(fdir);

cellsize = 0.125;
nodata_value = 0;

% Choose the appropriate definition:
xllcorner = min(lon);
yllcorner = min(lat);

georef = [xllcorner, yllcorner, cellsize, nodata_value];

% Do basin analysis. Includes snapping of gauges to flow map
snap = 1;
gage = zeros(m,nt); % placeholder
basin = basin_analysis_js(flipud(fdir), georef, gauge_list, gage, snap, lon, lat); % all gages
hold on 
mapshow(cropped_grwl(logical(id100)), 'SymbolSpec', symbspec, 'color', 'green')

figure,plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'flowacc');
hold on
plot(T2{:,1}, T2{:,2}, 'r.', 'markersize', 30)
mapshow(cropped_grwl(logical(id100)), 'SymbolSpec', symbspec, 'color', 'green')
title('Ohio basin with outlet location')

flow_vel = 1.4; % m/s
[HH, travel_time] = state_model(basin.mask, basin.mask_gage, basin.flow_length_gage, flow_vel, timestep);

cells_per_gage = sum(sum(HH,2),3);
k = size(HH, 3) - 1;

A = basin.mask(:);

[~, maxind] = max(basin.gage_area);

A(~isnan(A)) = travel_time{maxind};
B = flipud(reshape(A, nrows, ncols));

figure
plotraster(basin.lonv, basin.latv, ceil(B), 'Travel time (days)');

C = load('./ohio_data/ohio_nldas_nanfilled.mat');
basin_mask = flipud(basin.mask);
basin_mask(isnan(basin_mask)) = 0;
basin_mask_linear = basin_mask(:);
nldas_runoff_linear = reshape(C.nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_matrix = nldas_runoff_linear(logical(basin_mask_linear),:);
% save('./nldas_runoff_matrix.mat', 'nldas_runoff_matrix')
true_discharge = state_model_dumb(true_runoff', HH); % units are mm/day

figure,plot(true_discharge)

m = length(basin.gage_lat);

%% Sample the true discharge according to SWOT temporal sampling

% For each virtual gauge in basin.gage, find the nearest gauge location in
% T and use the associated overpass day to sample the discharge

orbit_day = zeros(m,1);
for i=1:m
    
    % Find nearest row in table
    a = [basin.gage_lon(i), basin.gage_lat(i)];
    b = [T2.lon, T2.lat];
    
    % Find a in b
    kk = dsearchn(b,a);
    
    orbit_day(i) = T2{kk,3};
    
end

true_discharge_w_swot_sampling = NaN(nt,m);
orbit_cycle = repmat(1:21,1,20);
orbit_cycle = orbit_cycle(1:nt);

for t=1:nt
    
    % get day of orbit cycle
    day_of_orbit_cycle = orbit_cycle(t);
    
    % get discharge for this time for locations matching the orbit day
    i1 = orbit_day == day_of_orbit_cycle;
    
    true_discharge_w_swot_sampling(t,i1) = true_discharge(t,i1);
        
end

figure,imagescnan(true_discharge_w_swot_sampling)
title('True discharge (mmd) (SWOT sampling)')
xlabel('Gage')
ylabel('Day')
colorbar
set(gca, 'fontsize', fs)

figure
plot(true_discharge(:,1),'linewidth', lw)
hold on
plot(true_discharge_w_swot_sampling(:,1), 'r.', 'markersize', ms)
title('Virtual gauge at basin outlet')
xlabel('Day of year')
ylabel('Discharge (mmd)')
legend('True discharge','SWOT measurement')
set(gca, 'fontsize', fs)

%% Save all data needed for ISR 

overpass_table = T2;

save('./ohio_data/swot_like_measurements_100m_no_error.mat', 'true_discharge',...
    'true_discharge_w_swot_sampling','overpass_table','basin','HH','k','orbit_cycle',...
    'flow_vel','travel_time','grwl100')
