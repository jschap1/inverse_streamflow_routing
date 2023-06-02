% Initial setup for Ohio basin ISR
%
% 6/2/2023 JRS
% 
% INPUTS:
% Flow directions file
% Start and end date (uses daily timestep)
%
% OUTPUTS:
% Identify measurement locations
% Basin analysis
% Measurement model

clean
addpath(genpath('./src/'))
addpath(genpath('/Users/jschap/Documents/Codes/VICMATLAB/VICMATLAB/vicmatlab'))

% Input start and end dates
tvec = datetime(2009,1,1):datetime(2009,12,31);

%% 1/8 flow directions and basin mask

fdfile = '/Volumes/HD3/ISR/05_Software/inverse_routing_mingpan/oh.dir';
load(fdfile); % oh
oh(oh==-1) = 0;

minlat = 34;
minlon = -90;
res = 0.125; % 1/8 degree resolution
R = makerefmat(minlon, minlat, res, res);
geotiffwrite('./basin_setup/ohio_flowdir.tif', oh, R);

fdfile = './basin_setup/ohio_flowdir.tif';
[fdir, R, lon, lat] = geotiffread2(fdfile);
fdir = flipud(fdir);

basin_mask = fdir;
basin_mask(basin_mask>0) = 1;
geotiffwrite('./basin_setup/ohio_basinmask.tif', basin_mask, R);

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

%% Identify all river channel pixels for ngage exps

figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'flowacc')

rivs = NaN(size(basin.flow_accum));
rivs(~isnan(basin.flow_accum)) = 0;
rivs(basin.flow_accum>1e10) = 1;
figure
plotraster(basin.lonv, basin.latv, flipud(rivs), 'River centerlines')

riv_pixel_coords = [basin.grid_lon(rivs == 1), basin.grid_lat(rivs == 1)];

m = length(riv_pixel_coords);
gauge_list = zeros(m,7);
for kk=1:m
    gauge_list(kk,:) = [1 1000+kk, riv_pixel_coords(kk,2) - 1/16, 360 + riv_pixel_coords(kk,1) - 1/16, 1000, 1000, 1000];
end

% make placeholder gage data
gage = zeros(m,nt);

% Do basin analysis. Includes snapping of gauges to flow map
snap = 0;
basin = basin_analysis_js(flipud(fdir), georef, gauge_list, gage, snap, lon, lat); % all gages

figure,plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'Ohio Basin with Gages');
hold on
plot(riv_pixel_coords(:,1)-1/16, riv_pixel_coords(:,2)-1/16, 'r.', 'markersize', 30)
title('Ohio basin with gage locations')

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

%% Plot number of downstream gages

% Plot map of number of downstream gauges for each cell
basin.ds_gauges = sum(basin.mask_gage,3);
figure
plotraster(basin.lonv, basin.latv, flipud(basin.ds_gauges), 'Number of downstream gauges')

%% Route NLDAS runoff to gage(s)

% Better to use routed runoff that is true to the NLDAS runoff and the flow
% direction file to avoid any inconsistencies introduced by USGS gauges

% Route NLDAS runoff to gauge locations
load('./ohio_nldas.mat')

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

%% Gauge discharge (not used right now...)

aa = load('/Volumes/HD3/ISR/05_Software/inverse_routing_mingpan/gage_03611500.txt');
figure, plot(tvec, aa(:,4)), title('USGS discharge')

% load discharge data from each gage
% gage_files = dir('/Volumes/HD3/ISR/05_Software/inverse_routing_mingpan/USGS/0*');
% nt = 365;
% m = 50;
% q = zeros(nt, m); % USGS gauge discharge (I think in units of cms)
% for i=1:m
%     tmp = load(fullfile('/Volumes/HD3/ISR/05_Software/inverse_routing_mingpan/USGS', gage_files(i).name));
%     q(:,i) = tmp(:,4);
% end

%% Save inputs

savedir = '/Volumes/HD3/ISR/02_Basins/Ohio/basin_setup';
save(fullfile(savedir, 'setup-392-gage.mat'),'basin','HH','flow_vel','k', 'gage', 'm', 'n', 'true_runoff', 'res', '-v7.3')


%% Generate all gage configurations of interest

% get flow accumulation associated with each river pixel
% sort by descending flow accumulation and take every other pixel
% to get every other pixel along the river channel
%
% It is not perfect, but not terrible, either

[basin.grid_lon(rivs == 1), basin.grid_lat(rivs == 1)];
riv_pixel_facc = basin.flow_accum(rivs == 1);
[B,I] = sort(riv_pixel_facc, 'descend');
dlmwrite('./riv_pixel_facc.txt', riv_pixel_facc)

dlmwrite('./river_pixel_coordinates.txt', riv_pixel_coords)

rpc_sort = riv_pixel_coords(I,:);

% figure

basin_orig = basin;

for inc = 2:100
%     subplot(5,5,inc)
    grid_length = 12.5; % km, approx
    mi = 1:inc:392;
%     plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'Ohio Basin with Gages');
%     hold on
%     plot(rpc_sort(mi,1)-1/16, rpc_sort(mi,2)-1/16, 'r.', 'markersize', 30)
%     title(['Approx. spacing: ' num2str(inc*grid_length)])
    
    qq = I(mi);
    basin.gage_lat = basin_orig.gage_lat(qq);
    basin.gage_lon = basin_orig.gage_lon(qq);
    basin.gage_i = basin_orig.gage_i(qq);
    basin.gage_j = basin_orig.gage_j(qq);
    basin.gage_area = basin_orig.gage_area(qq);
    basin.gage_accum = basin_orig.gage_accum(qq);
    basin.mask_gage = basin_orig.mask_gage(:,:,qq);
    basin.flow_length_gage = basin_orig.flow_length_gage(:,:,qq);
    basin.ia = basin_orig.ia(qq);
    basin.gage = basin_orig.gage(qq,:);
    basin.ds_gauges = sum(basin.mask_gage,3);
    basin.spacing = inc*grid_length;
    
    flow_vel = 1.4; % m/s
    [HH, travel_time] = state_model(basin.mask, basin.mask_gage, basin.flow_length_gage, flow_vel, timestep);

    cells_per_gage = sum(sum(HH,2),3);
    basin.cells_per_gage = cells_per_gage;
    k = size(HH, 3) - 1;

    discharge = state_model_dumb(nldas_runoff_matrix', HH); % units are mm/day
    gage = discharge;
    true_runoff = nldas_runoff_matrix;

    save(fullfile(savedir, ['setup-' num2str(length(basin.gage_lon)) '-gage.mat']),...
        'basin','HH','flow_vel','k', 'gage', 'm', 'n', 'true_runoff', 'res', '-v7.3')

    disp(['Generated inputs for ' num2str(inc*grid_length) ' km spacing'])
end

figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'Ohio Basin with Gages');
hold on
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', 30)
title('Ohio basin with gage locations')
