% Initial setup for Hetch Hetchy basin ISR
%
% 6/14/2023 JRS
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
addpath(genpath('./src/'))

%% Flow directions

fdfile = '/hdd/ISR/02_Basins/HH/inverse_routing_dev/flowdir_vic.tif';
[fdir, R, lon, lat] = geotiffread2(fdfile);

fdir = double(fdir(:,2:end));
fdir(isnan(fdir)) = 0;
fdir(3,1) = 9;

figure,plotraster(lon, lat, fdir, 'flow directions');

%% Load true runoff

load('/hdd/ISR/02_Basins/HH2/Data/unsorted/truth.mat');
nt = length(truth.t);

%% Basin analysis

timestep = 3600; % hourly timestep
[nrows, ncols] = size(fdir);
n = sum(fdir(:) ~= 0);

cellsize = 1/16;
nodata_value = 0;

% Choose the appropriate definition:
xllcorner = min(lon);
yllcorner = min(lat);

georef = [xllcorner, yllcorner, cellsize, nodata_value];

gage_list = load('./hh_data/station_list_5.txt');
m = size(gage_list,1);

% Do basin analysis. Includes snapping of gauges to flow map
snap = 0;
gage = zeros(m, nt); % placeholder
basin = basin_analysis_js(flipud(fdir), georef, gage_list, gage, snap, lon, lat); % all gages

basin.gage_area(1)/1e6
basin.fdir = fdir;

%% Visualize in Flow Direction Toolkit

% flow direction
% /hdd/ISR/02_Basins/HH/inverse_routing_dev/flowdir_vic.tif

% basin boundary
% /hdd/ISR/02_Basins/HH/Forward_Modeling/tuosub_bb.shp

% gages
% /hdd/ISR/02_Basins/HH2/Data/gage_loc_5.txt

% river
% /hdd/ISR/02_Basins/HH/Forward_Modeling/Delineation/tuorivs_crop.shp

% Remove rivers that fall outside of tuosub_bb
% ogr2ogr -t_srs epsg:4326 -s_srs epsg:4326 -clipsrc /hdd/ISR/02_Basins/HH/Forward_Modeling/tuosub_bb.shp /hdd/ISR/02_Basins/HH/Forward_Modeling/Delineation/tuorivs_crop.shp /hdd/ISR/02_Basins/HH/Forward_Modeling/Delineation/usrivs_crop.shp

%% Create state transition model

flow_vel = 2; % m/s
[HH, travel_time] = state_model(basin.mask, basin.mask_gage, basin.flow_length_gage, flow_vel, timestep);
cells_per_gage = sum(sum(HH,2),3);
k = size(HH, 3) - 1;

A = basin.mask(:);
[~, maxind] = max(basin.gage_area);
A(~isnan(A)) = travel_time{maxind}; % off by 1 cell
B = flipud(reshape(A, nrows, ncols));
figure
plotraster(basin.lonv, basin.latv, ceil(B), 'Travel time (hours)');

%% Generate true discharge at gauge locations

truth.discharge_mm = state_model_dumb(truth.total_runoff', HH);

%% Save ISR inputs

save('./hh_data/isr_setup_usds_gage_5.mat', 'basin', 'truth', 'flow_vel', 'k', 'HH', 'fdir')
