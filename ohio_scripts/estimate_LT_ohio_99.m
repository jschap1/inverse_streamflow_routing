% Estimating L, T for Ohio basin
%
% 6/12/2023 JRS
% Uses our knowledge of the true runoff to estimate L,T
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
cd /home/jschap/Documents/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

load('./ohio_data/ohio_tmpa_nanfilled.mat')
load('./ohio_data/ohio_nldas_nanfilled.mat')
load('./ohio_data/swot_like_measurements_100m_no_error_revised.mat')

% add distmat
A = load('./ohio_data/setup-1-gage.mat');
basin.distmat = A.basin.distmat;

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

%% Calculate true runoff error and export to R (easier to analyze there)

figure
subplot(1,2,1)
plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,1), 'TMPA (day 1)')
subplot(1,2,2)
plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,1), 'NLDAS (day 1)')

runoff_error = nldas.runoff - tmpa.runoff; % truth - estimate

% what is the runoff error temporal decorrelation length?
% what is the runoff error spatial decorrelation length?

% reshape to 2D matrix for writing to file
[nr,nc,nt] = size(runoff_error);
runoff_error_2d = reshape(runoff_error, nr*nc, nt);
dlmwrite('./ohio_data/nldas_tmpa_error.txt', runoff_error_2d, 'delimiter', '\t')
dlmwrite('./ohio_data/distmat.txt', basin.distmat, 'delimiter', '\t')

for tt=1:365
    [x,y,z] = grid2xyz(basin.lonv', basin.latv', runoff_error(:,:,tt));
    dlmwrite(['./ohio_data/runoff_error_xyz/xyz' num2str(tt) '.txt'], [x,y,z], 'delimiter', '\t')
end