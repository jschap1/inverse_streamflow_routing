% Initial setup for Ohio basin ISR
%
% 6/5/2023 JRS
% 
% INPUTS:
% Flow directions file with 9 at outlet
% NLDAS runoff data
% Start and end date (uses daily timestep)
%
% OUTPUTS:
% Identify measurement locations
% Basin analysis
% Measurement model

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

% Input start and end dates
tvec = datetime(2009,1,1):datetime(2009,12,31);

%% 1/8 flow directions and basin mask

fdfile = './ohio_data/oh.dir';
load(fdfile); % oh
oh(oh==-1) = 0;

minlat = 34;
minlon = -90;
res = 0.125; % 1/8 degree resolution
R = makerefmat(minlon, minlat, res, res);
geotiffwrite('./ohio_data/ohio_flowdir.tif', oh, R);

fdfile = './ohio_data/ohio_flowdir.tif';
[fdir, R, lon, lat] = geotiffread2(fdfile);
fdir = flipud(fdir);

basin_mask = fdir;
basin_mask(basin_mask>0) = 1;
geotiffwrite('./ohio_data/ohio_basinmask.tif', basin_mask, R);

basin_mask_nan = basin_mask;
basin_mask_nan(basin_mask_nan==0) = NaN;
n = sum(basin_mask(:)); % number of grid cells

figure,plotraster(lon, lat, fdir, 'flow directions');
figure,plotraster(lon, lat, basin_mask, 'basin mask');

%% Make gage list

% for a single gauge at the outlet
outlet_lon = -89.25;
outlet_lat = 37;
gauge_list = [1 11111111, outlet_lat, 360 + outlet_lon, 1000, 1000, 1000];

m = size(gauge_list,1);

gagecoords = [gauge_list(:,4) - 360, gauge_list(:,3)];

figure,plotraster(lon, lat, basin_mask, 'basin mask');
hold on
plot(gagecoords(:,1), gagecoords(:,2), 'r.', 'markersize', 30)
title('Ohio basin with outlet location')

% make placeholder gage data
nt = 365;
gage = zeros(m,nt);

%% Basin prep

timestep = 86400; % s
[nrows, ncols] = size(fdir);

cellsize = 0.125;
nodata_value = 0;

% Choose the appropriate definition:
xllcorner = min(lon);
yllcorner = min(lat);

georef = [xllcorner, yllcorner, cellsize, nodata_value];

% Do basin analysis. Includes snapping of gauges to flow map
snap = 0;
basin = basin_analysis_js(flipud(fdir), georef, gauge_list, gage, snap, lon, lat); % all gages

figure,plotraster(basin.lonv, basin.latv, fdir, 'flowdir');
hold on
plot(gagecoords(:,1), gagecoords(:,2), 'r.', 'markersize', 30)
title('Ohio basin with outlet location')

%% Create routing (state transition) model

flow_vel = 1.4; % m/s
[HH, travel_time] = state_model(basin.mask, basin.mask_gage, basin.flow_length_gage, flow_vel, timestep);

cells_per_gage = sum(sum(HH,2),3);
k = size(HH, 3) - 1;

A = basin.mask(:);

[~, maxind] = max(basin.gage_area);

A(~isnan(A)) = travel_time{maxind}; % off by 1 cell
B = flipud(reshape(A, nrows, ncols));

figure
plotraster(basin.lonv, basin.latv, ceil(B), 'Travel time (days)');

%% Plot number of downstream gages for each cell

basin.ds_gauges = sum(basin.mask_gage,3);
figure
plotraster(basin.lonv, basin.latv, flipud(basin.ds_gauges), 'Number of downstream gauges')

%% Route NLDAS runoff to gage(s)

% Better to use routed runoff that is true to the NLDAS runoff and the flow
% direction file to avoid any inconsistencies introduced by USGS gauges

% Route NLDAS runoff to gauge locations
load('./ohio_data/ohio_nldas.mat')

% Fill missing values in the NLDAS runoff data
aa = nldas.runoff(:,:,1);
aa(isnan(aa)) = 0;
aa(basin_mask==0) = NaN;
[i1,i2] = find(aa==0); % index of missing values
nldas.runoff(i1,i2,:) = mean([nldas.runoff(i1+1,i2,:), nldas.runoff(i1-1,i2,:), ...
    nldas.runoff(i1,i2+1,:), nldas.runoff(i1,i2-1,:)]);

save('./ohio_data/ohio_nldas_nanfilled.mat', 'nldas')
    
% Repeat for TMPA runoff data
load('./ohio_data/ohio_tmpa.mat')
aa = tmpa.runoff(:,:,1);
aa(isnan(aa)) = 0;
aa(basin_mask==0) = NaN;
[i1,i2] = find(aa==0); % index of missing values
tmpa.runoff(i1,i2,:) = mean([tmpa.runoff(i1+1,i2,:), tmpa.runoff(i1-1,i2,:), ...
    tmpa.runoff(i1,i2+1,:), tmpa.runoff(i1,i2-1,:)]);
save('./ohio_data/ohio_tmpa_nanfilled.mat', 'tmpa')

figure, subplot(1,2,1)
plotraster(lon,lat, nldas.runoff(:,:,1), 'nldas runoff')
subplot(1,2,2)
plotraster(lon,lat, nldas.runoff(:,:,1).*basin_mask, 'masked')

basin_mask_linear = basin_mask(:);
nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_matrix = nldas_runoff_linear(logical(basin_mask_linear),:);
discharge = state_model_dumb(nldas_runoff_matrix', HH); % units are mm/day

figure
plot(tvec, discharge), title('NLDAS discharge at outlet')
xlabel('Time'), ylabel('Discharge (mm/day)')

gage = discharge;
true_runoff = nldas_runoff_matrix;

%% Calculate distance matrix (can take some time)

aaa = load('/Volumes/HD3/ISR/02_Basins/Ohio/basin_setup/setup-1-gage.mat');
basin.distmat = aaa.basin.distmat;

% basinmask = basin.mask;
% basinmask(isnan(basinmask)) = 0;
% basinmask = logical(basinmask);
% basin.distmat = calc_distance_matrix(basin.grid_lon(basinmask), basin.grid_lat(basinmask));

%% Save inputs

save('./ohio_data/setup-1-gage.mat','basin','HH','flow_vel','k', 'gage', 'm', 'n', 'true_runoff', 'res', '-v7.3')
