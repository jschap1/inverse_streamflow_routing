% Initial setup for Allegheny basin ISR
%
% 7/17/2023 JRS
% 
% INPUTS:
% Flow directions file with 9 at outlet
% NLDAS runoff data
% Start and end date (uses daily timestep)
% Uses Ohio basin setup and subsets it, basically
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

% Load mask for cropping Ohio to Alleg subbasin
mask = geotiffread('./allegheny_data/alleg_mask.tif');
mask = flipud(mask);
[keeprows, keepcols] = find(mask==1);
keeprows = unique(keeprows);
keepcols = unique(keepcols);
alleg_mask = mask(keeprows,keepcols);
mask(mask==0) = NaN;

figure
imagesc(alleg_mask)

%% 1/8 flow directions and basin mask

fdfile = './ohio_data/oh.dir';
load(fdfile); % oh
oh(oh==-1) = 0;

minlat = 40;
minlon = -80.38;
res = 0.125; % 1/8 degree resolution
R = makerefmat(minlon, minlat, res, res);

fdfile = './ohio_data/ohio_flowdir.tif';
[fdir, ~, lon, lat] = geotiffread2(fdfile);
fdir = flipud(fdir);

figure
subplot(2,1,1)
imagesc(mask), title('mask')
subplot(2,1,2)
imagesc(fdir), title('fdir')

% Crop fdir to Alleg
fdir = fdir.*(mask);
fdir = fdir(keeprows, keepcols);

geotiffwrite('./allegheny_data/alleg_basinmask.tif', alleg_mask, R);

n = nansum(alleg_mask(:)); % number of grid cells

figure,plotraster(lon, lat, fdir, 'flow directions');
figure,plotraster(lon, lat, alleg_mask, 'basin mask');

alleg_mask(alleg_mask==0) = NaN;
fdir = fdir.*alleg_mask;

figure,plotraster(lon, lat, fdir, 'flow directions');
figure,plotraster(lon, lat, alleg_mask, 'basin mask');

geotiffwrite('./allegheny_data/allegheny_fd_cropped.tif', fdir, R)
geotiffwrite('./allegheny_data/allegheny_mask_cropped.tif', alleg_mask, R)

[~,~,lon, lat] = geotiffread2('./allegheny_data/allegheny_mask_cropped.tif');

fdir = geotiffread2('./allegheny_data/allegheny_fd_cropped.tif');

%% Some gauge locations of interest

% Allegheny River at Natrona, PA 11400 sq. miles. (A = 29525 km2)
[-79.7186, 40.6153] % NAD27
[-79.7186000, 40.6153000] % WGS84

% Allegheny River near Rimer, PA 8389 sq. miles (A = 21727 km2)
[-79.5472, 40.9583] % NAD27
[-79.5472000, 40.9583000] % WGS84

% Mahoning Creek at Mahoning Creek Dam, PA 344 sq. miles (A = 891 km2)
[-79.2914, 40.9275] % NAD27
[-79.2914000, 40.9275000] % WGS84

%% Make gage list

A = load('./ohio_data/swot_like_measurements_100m_no_error_revised.mat');
g_lon = A.basin.gage_lon;
g_lat = A.basin.gage_lat;

% g_lon = [-79.7186, -79]; % orig
% g_lat = [40.6153, 41.75]; % orig
g_lon = [-79.75]; % outlet
g_lat = [40.5]; % outlet
% g_lon = [-79.75, -79]; % revised
% g_lat = [40.5, 41.75]; % revised
% g_lon = [-79.7186, -79.5472, -79.2914];
% g_lat = [40.6153, 40.9583, 40.9275];
for i=1:length(g_lon)
    gauge_list(i,:) = [1, i, g_lat(i), 360 + g_lon(i), 1000, 1000, 1000];
end

% gagecoords(1,:) = [];
% gagecoords(233,:) = [-80.2550-res, 41.75];

gagecoords = [gauge_list(:,4) - 360, gauge_list(:,3)];
% gagecoords = [gauge_list(:,4) - 360-res/2, gauge_list(:,3)-res/2];

figure,plotraster(lon, lat, alleg_mask, 'basin mask');
hold on
plot(gagecoords(:,1), gagecoords(:,2), 'r.', 'markersize', 30)
title('Allegheny basin with gages')

figure,plotraster(lon, lat, fdir, 'basin mask');
hold on
plot(gagecoords(:,1), gagecoords(:,2), 'r.', 'markersize', 30)
title('Allegheny flow directions basin with outlet location')

% Just keep the gauges that are in the Allegheny basin
rm1=gagecoords(:,1)<min(lon);
rm2=gagecoords(:,2)<min(lat);
rm3=gagecoords(:,1)>max(lon);
rm4=gagecoords(:,2)>max(lat);
gagecoords(rm1|rm2|rm3|rm4,:)=[];
gauge_list(rm1|rm2|rm3|rm4,:)=[];
m = size(gagecoords,1);

% make placeholder gage data
nt = 365;
gage = zeros(m,nt);

% set fdir=9 at outlet
fdir(5,6) = 9; 

% estimated basin size: 
% 1) 33750 km2
% 2) 20870 km2
% 3) 5233 km2

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

% gagecoords(1,:) = [];
% gagecoords(233,:) = [-80.2550-res, 41.75];

% bb=basin.mask;
% bb(isnan(bb))=0;
% bb=logical(bb);
% lonvec = basin.grid_lon(bb);
% latvec = basin.grid_lat(bb);
% % generate 10 random gauges
% rng(704753262, 'twister')
% gagei=1:n;
% % gagei = randperm(233,9);
% g_lon = [lonvec(gagei)];
% g_lat=[latvec(gagei)];
% % g_lon = [g_lon;lonvec(gagei)];
% g_lat=[g_lat;latvec(gagei)];

%% Create routing (state transition) model

flow_vel = 1.4; % m/s
[HH, travel_time] = state_model(basin.mask, basin.mask_gage, basin.flow_length_gage, flow_vel, timestep);

cells_per_gage = sum(sum(HH,2),3);
k = size(HH, 3) - 1;

A = basin.mask(:);

[~, maxind] = max(basin.gage_area);

[nrows, ncols] = size(basin.mask);
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

% Route NLDAS/TMPA runoff to gauge locations
load('./ohio_data/ohio_nldas_nanfilled.mat')
load('./ohio_data/ohio_tmpa_nanfilled.mat')

alleg_mask(alleg_mask==0) = NaN;

figure, imagesc(mask)
figure, imagesc(nldas.runoff(:,:,1));

nldas.runoff = nldas.runoff(keeprows,keepcols,:);
nldas.runoff = nldas.runoff.*alleg_mask;

nldas.prec = nldas.prec(keeprows,keepcols,:);
nldas.prec = nldas.prec.*alleg_mask;

tmpa.runoff = tmpa.runoff(keeprows,keepcols,:);
tmpa.runoff = tmpa.runoff.*alleg_mask;

tmpa.prec = tmpa.prec(keeprows,keepcols,:);
tmpa.prec = tmpa.prec.*alleg_mask;

% save('./allegheny_data/alleg_nldas_nanfilled.mat', 'nldas')
% save('./allegheny_data/alleg_tmpa_nanfilled.mat', 'tmpa')

figure,
plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,222), 'tmpa')

figure
subplot(1,2,1)
plotraster(lon,lat, nldas.runoff(:,:,1), 'nldas runoff')
subplot(1,2,2)
plotraster(lon,lat, tmpa.runoff(:,:,1), 'tmpa runoff')

figure
subplot(1,2,1)
plotraster(lon,lat, nldas.prec(:,:,1), 'nldas prec')
subplot(1,2,2)
plotraster(lon,lat, tmpa.prec(:,:,1), 'tmpa prec')

basin.mask(isnan(basin.mask)) = 0;
basin.mask = flipud(basin.mask);
basin_mask_linear = basin.mask(:);

nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_matrix = nldas_runoff_linear(logical(basin_mask_linear),:);
discharge = state_model_dumb(nldas_runoff_matrix', HH); % units are mm/day

figure
plot(tvec, discharge), title('NLDAS discharge')
xlabel('Time'), ylabel('Discharge (mm/day)')

gage = discharge;
true_runoff = nldas_runoff_matrix;

%% Calculate distance matrix (can take some time)

% number of cells in each subbasin
n1 = sum(nansum(basin.mask_gage(:,:,1)));
n2 = sum(nansum(basin.mask_gage(:,:,2)));

% aaa = load('/Volumes/HD3/ISR/02_Basins/Ohio/basin_setup/setup-1-gage.mat');
% basin.distmat = aaa.basin.distmat;

basinmask = flipud(basin.mask);
basinmask(isnan(basinmask)) = 0;
basinmask = logical(basinmask);
basin.distmat = calc_distance_matrix(basin.grid_lon(basinmask), basin.grid_lat(basinmask));

figure
imagescnan(basin.distmat), title('distance matrix (km)')

%% Save inputs

A = load('/Volumes/HD3/ISR/02_Basins/Ohio/basin_setup/setup-1-gage.mat');

% map of Allegheny in context of Ohio
figure
plotraster(lon, lat, alleg_mask, 'Study area')
hold on
figure
plotraster(A.basin.lonv, A.basin.latv, flipud(A.basin.mask), 'Study area')
hold on
imagesc(lon, lat, alleg_mask)

tmpa_runoff = tmpa.runoff;
nldas_runoff = nldas.runoff;

save('./allegheny_data/setup-swot-gage.mat','basin','HH','flow_vel','k', ...
    'gage', 'm', 'n', 'true_runoff', 'res', 'nldas_runoff','tmpa_runoff','-v7.3')

%% Stats


n = sum(nansum(basin.mask_gage(:,:,1)));
sqrt_basin_area = sqrt(sum(nansum(basin.gage_area(1)))/1e6); % km2
time_of_conc = max(nanmax(basin.flow_length_gage(:,:,1)))/flow_vel/86400; % days
