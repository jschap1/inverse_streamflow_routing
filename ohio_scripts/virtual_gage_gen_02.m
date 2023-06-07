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
symbspec = makesymbolspec('Line',{'Default','Color','#7E2F8E'});

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


%% Get river widths

% Make a bounding box shapefile for cropping purposes

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

%% Create gauge list and route true runoff to each SWOT overpass location

T3 = dlmread('./ohio_data/overpass_table.txt', '\t', 1, 0);
% T3 = readmatrix('./ohio_data/overpass_table.txt');
T3 = table(T3(:,1), T3(:,2), T3(:,3));
T3.Properties.VariableNames{1} = 'lon';
T3.Properties.VariableNames{2} = 'lat';
T3.Properties.VariableNames{3} = 'orbit_day';

virtual_gage_locs = unique([T3.lon, T3.lat], 'rows');
m = size(virtual_gage_locs,1);

nt = 365;
gauge_list = [ones(m,1), 1000+(1:m)', virtual_gage_locs(:,2), virtual_gage_locs(:,1) + 360, 1000*ones(m,3)];

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

figure,plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'flowacc');
hold on
plot( virtual_gage_locs(:,1),  virtual_gage_locs(:,2), 'r.', 'markersize', 30)
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

m = length(basin.gage_lat); % There are 244 virtual gauges for 100m 
% some of these gauges have multiple observations

%% Sample the true discharge according to SWOT temporal sampling

% For each virtual gauge in basin.gage, find the nearest gauge location in
% T and use the associated overpass day to sample the discharge

orbit_days = cell(m,1);
isr_res = 1/8; 
for i=1:m
    
    % Find nearest row in table
    a = [basin.gage_lon(i), basin.gage_lat(i)];
    b = [T3.lon, T3.lat];
    
    % Find all observations for a gage
    kk = find(abs(T3.lon - basin.gage_lon(i))<=isr_res & abs(T3.lat - basin.gage_lat(i))<=isr_res);
    
    % Find a in b
%     kk = dsearchn(b,a);
    
    orbit_days{i} = unique(T3{kk,3});
    
end

find(orbit_days==0)
% gages have between 0-4 observations per cycle

%% Get discharge for this time for locations matching the orbit day

true_discharge_w_swot_sampling = NaN(nt,m);
orbit_cycle = repmat(1:21,1,20);
orbit_cycle = orbit_cycle(1:nt);

for t=1:nt
    
    % get day of orbit cycle
    day_of_orbit_cycle = orbit_cycle(t);
    
    % get discharge for this time for locations matching the orbit day
    i1 = [];
    for mm=1:m
        if any(orbit_days{mm} == day_of_orbit_cycle)
            i1 = [i1; mm];
        end
    end
    
    % index of which gages are observed on this day
%     i1 = orbit_days == day_of_orbit_cycle;
    
    true_discharge_w_swot_sampling(t,i1) = true_discharge(t,i1);
        
end

% Remove gauges with no observations:
i_remove = find(num_overpasses==0);
true_discharge(:,i_remove) = [];
true_discharge_w_swot_sampling(:,i_remove) = [];
orbit_days(i_remove) = [];

%% Check virtual measurements

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

overpass_table = T3;

save('./ohio_data/swot_like_measurements_100m_no_error_revised.mat', 'true_discharge',...
    'true_discharge_w_swot_sampling','overpass_table','basin','HH','k','orbit_cycle',...
    'flow_vel','travel_time')

%% Make map of virtual gage locations for paper

m = length(orbit_days);
for mm=1:m
    num_overpasses(mm) = length(orbit_days{mm});
end

grwl = shaperead('./ohio_data/grwl_ohio.shp');
swot = shaperead('./ohio_data/swot_ohio.shp');

figure
plotraster(basin.lonv, basin.latv, flipud(basin.mask), 'Virtual Gages (>100m)')
colorbar off
hold on
for i=1:length(grwl)
    plot(grwl(i).X, grwl(i).Y, 'k') % grwl centerlines
end
% for i=1:length(swot)
%     plot(swot(i).X, swot(i).Y, 'blue') % swot overpasses
% end
numGroups = length(unique(num_overpasses));
clr = [254,240,217;
253, 204, 13;
227, 74, 51; 
179, 0, 0]/255;
siz = [1,1,1,1,1]*30;
h = gscatter(basin.gage_lon, basin.gage_lat, num_overpasses', clr, '.', siz, 'off'); % color code by nobs/cycle
legend(h, 'n_{obs} = 1', 'n_{obs} = 2', 'n_{obs} = 3', 'n_{obs} = 4')
axis('equal')



