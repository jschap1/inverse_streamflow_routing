%% Running Y20 ISR
%
% 10/5/2023 JRS
% Implementing domain Y20 ISR for Allegheny basin

clear, clc, close all
cd /home/jschap/Dropbox/ISR/inverse_streamflow_routing/
% cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))
lw = 2;

A=load('./allegheny_data/setup/setup-2-gage.mat');
load('./allegheny_data/setup/setup-swot-gage.mat');

B = load('./allegheny_data/setup/setup-233-gage.mat');
basin = B.basin;
gage = B.gage;
HH = B.HH;

basin.distmat = A.basin.distmat;
clearvars A B
load('./allegheny_data/setup/alleg_nldas_nanfilled.mat')
load('./allegheny_data/setup/alleg_tmpa_nanfilled.mat')

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);
true_discharge = gage;
[nt,m] = size(true_discharge);

%% Cut down to a shorter time period

i1=60; i2=179;
tv = tv(i1:i2);
gage = gage(i1:i2,:);
nldas_runoff = nldas_runoff(:,:,i1:i2);
true_runoff = true_runoff(:,i1:i2); 
[nt,m] = size(gage);
nldas.runoff  = nldas.runoff(:, :, i1:i2);
tmpa.runoff  = tmpa.runoff(:, :, i1:i2);
true_discharge = true_discharge(i1:i2,:);
nt=length(tv);
basin.true_runoff = true_runoff';

%% Reshape true runoff from map to linear

% basin.mask = flipud(basin.mask);
figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be right side up

basin_mask_linear = basin.mask(:);
basin_mask_linear(isnan(basin_mask_linear)) = 0;

tmpa_runoff_linear = reshape(tmpa.runoff, length(basin_mask_linear), nt);
tmpa_runoff_prior = tmpa_runoff_linear(logical(basin_mask_linear),:)';
tmpa_runoff_prior(isnan(tmpa_runoff_prior)) = 0;

nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_true = nldas_runoff_linear(logical(basin_mask_linear),:)';
nldas_runoff_true(isnan(nldas_runoff_true)) = 0;

%% Generate discharge "measurements" (Y20)

rng(704753262,'twister')

% Y20: additive Gaussian error with mu, sigma
swot = 0;
freq = 10;
percent_error = 5;

mu1 = 0; % mean of additive error
% sigma1 = 0.3*gage; % stddev of additive error (30% of truth)
add_err = zeros(nt, m);
for tt=1:nt
    for mm=1:m
        sigma1 = (percent_error/100)*gage(tt,mm);
        add_err(tt,mm) = mu1 + sigma1*randn(1,1);
    end
    disp(tt)
end
gage_w_error = gage + add_err;

if swot  
   tmp = gage_w_error;
   gage_w_error = NaN(size(tmp));
   gage_w_error(1:freq:end,:) = tmp(1:freq:end,:);
end

%% Set aside validation gauges 

rng(704753262, 'twister');
nval = 7; 
val_ind = 1 - 1 + randperm(14-1+1,nval); % setting aside nval gauges for validation
cal_ind = find(~ismember(1:14, val_ind));

% two gauges
cal_ind = [29];
val_ind = [146];

% outlet only
%cal_ind = 2;
%val_ind = [1,3:14];

% all gauges
% cal_ind = 1:14;
% val_ind = [];

ms = 10;
res = 1/8;
figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'Ohio basin')
hold on
plot(basin.gage_lon(cal_ind)-res/2, basin.gage_lat(cal_ind)-res/2, 'r.', 'markersize', ms)
plot(basin.gage_lon(val_ind)-res/2, basin.gage_lat(val_ind)-res/2, 'green.', 'markersize', ms)
legend('', 'Calibration gauges','Validation gauges')

%% Plot discharge "measurements" (Y20)

ms = 30;
fs = 16;
figure
plot(tv, true_discharge(:,cal_ind(1)), 'linewidth', lw)
hold on 
plot(tv, gage_w_error(:,cal_ind(1)), '.', 'markersize', ms)
legend('Truth','Measurements')
xlabel('time')
ylabel('Q (mm/day)')
grid on
title('Discharge at outlet')
set(gca, 'fontsize', fs)

%% Load runoff prior (additive, ensemble of priors)

% /m0aflatL40T5_AG/runoff_errors_for_variogram_AGflat.mat
a = load('./allegheny_data/errors/m0a1L40T5/Ens_errors_AG.mat');
runoff_errors = a.runoff_errors;
[n,nt,M] = size(runoff_errors);
clearvars a
% Generate prior
prior_runoff_ens = zeros(n, nt, M);
for mm=1:M
    prior_runoff_ens(:,:,mm) = true_runoff + runoff_errors(:,:,mm);
end

%% Do ISR on each prior realization, using Y20 method

total_mean = mean(true_runoff(:));
% corrfact = mean(true_runoff(:))/mean(tmpa_runoff_prior(:));
prior_runoff_ens = tmpa_runoff_prior';
% prior_runoff_ens = corrfact*tmpa_runoff_prior'; % medium value
% prior_runoff_ens = 0.25*corrfact*tmpa_runoff_prior'; % low value
% prior_runoff_ens = 4*corrfact*tmpa_runoff_prior'; % high value
% M = size(prior_runoff_ens,3);

% uniform priors
% prior_runoff_ens = total_mean*ones(n,nt)./4; % low
% prior_runoff_ens = total_mean*ones(n,nt); % medium
% prior_runoff_ens = 4*total_mean*ones(n,nt); % high

%M = 2000; % total number of replicates
M = 1;
addgauss = 1;
multlog = 0;
if multlog
    a = load('./allegheny_data/errors/m1a1L40T5_LM/simulated_runoff_errors_2000.mat');
    runoff_errors = reshape(a.alldata, 20*21, nt, M);
    basinmask = basin.mask; % no nan values
    basinmask(isnan(basin.mask)) = 0;
    basinmask = logical(basinmask);
    runoff_errors = runoff_errors(basinmask, :, :);
end

%M=500;
prior_runoff = zeros(nt,n,M);
post_runoff = zeros(nt,n,M); % stores of M each runoff priors
prior_stddev = zeros(nt,n,M);
post_stddev = zeros(nt,n,M);
post_stddevQ = zeros(nt,m,M);
prior_stddevQ = zeros(nt,m,M);

tic

for mm=1:M

% Loading the runoff prior   

% mm=1;
if addgauss
    % add/gauss
    runoff_prior = prior_runoff_ens(:,:,mm)'; % may need to transpose
elseif multlog
    % mm=1; % just using one replicate
    % mult/log
    runoff_prior = nldas_runoff_true'.*runoff_errors(:,:,mm);
    runoff_prior = runoff_prior';
    prior_runoff(:,:,mm) = runoff_prior;
end

% The Actual ISR Part

s = 2*(k+1);
optionsfile = './allegheny_data/setup/options_alleg.txt';
rho_thres = 0;
L = 40; % km
T = 5; % days
alpha1 = 0.1;

[post_runoff(:,:,mm), post_stddev(:,:,mm), prior_stddev(:,:,mm)] = ISR_Y20(runoff_prior, HH(cal_ind,:,:), ...
    gage_w_error(:,cal_ind), s, basin, optionsfile, L, T, rho_thres, alpha1);

disp(['Finished run: ' num2str(mm)])

end % end loop over prior realizations

toc

if multlog
    prior_runoff_ens = prior_runoff;
end

%% Save results
save('./allegheny_data/results/E6_Y20AG_m0a1L40T5_daily_5_tmpa_alpha0.1.mat', 'runoff_prior', 'post_runoff', ...
    'gage_w_error','post_stddev', 'prior_stddev', 'post_stddevQ','prior_stddevQ',...
    'tmpa_runoff_prior', 'nldas_runoff_true', 'prior_runoff_ens', 'cal_ind', 'val_ind');

% may need to show a very optimistic measureemnt case to show how mult
% prior causes issues.. (11/18/23)

% save('./allegheny_data/results/allgage/Y20AGflat_m0L40T5_swot_15_1-M.mat', 'runoff_prior', 'post_runoff', ...
%     'gage_w_error','post_stddev', 'prior_stddev');
% save('./allegheny_data/errors/Y20AGflat_m0L40T5_prior_ensemble.mat', 'prior_runoff_ens')
