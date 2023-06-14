%% Comparing Y20 ISR and ensemble ISR
%
% 6/5/2023 JRS

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

load('./ohio_data/ohio_tmpa_nanfilled.mat')
load('./ohio_data/ohio_nldas_nanfilled.mat')
load('./ohio_data/swot_like_measurements_100m_no_error_revised.mat')

% add distmat
A = load('./ohio_data/setup-1-gage.mat');
basin.distmat = A.basin.distmat;

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

% load('./basin_setup/basin.mat')

%% Crop to a shorter domain

% tmax = 50;
% tv = tv(1:tmax);
% true_discharge = true_discharge(1:tmax,:);
% true_discharge_w_swot_sampling = true_discharge_w_swot_sampling(1:tmax,:);
% tmpa.runoff = tmpa.runoff(:,:,1:tmax);
% nldas.runoff = nldas.runoff(:,:,1:tmax);

[nt,m] = size(true_discharge);

%% Plot tmpa prior vs. nldas (true) runoff (Y20 Figure 6)

% TMPA prior, NLDAS truth

figure
t=43:47;
for i=1:5
    subplot(2,5,i)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['TMPA runoff (day ' num2str(t(i)) ')'])
    caxis([0,6])
    subplot(2,5,5+i)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['NLDAS runoff (day ' num2str(t(i)) ')'])
    caxis([0,6])
end
colormap cool

%% Assemble runoff prior (Y20)

basin.mask = flipud(basin.mask);
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

% gage = true_discharge; % should corrupt with error
gage = true_discharge_w_swot_sampling; % should corrupt with error

% Y20: additive Gaussian error with mu, sigma

mu1 = 0; % mean of additive error
% sigma1 = 0.3*gage; % stddev of additive error (30% of truth)
for tt=1:nt
    for mm=1:m
        sigma1 = 0.15*gage(tt,mm);
        add_err(tt,mm) = mu1 + sigma1*randn(1,1);
    end
    disp(tt)
end
gage_w_error = gage + add_err;
% add_err = mu1 + sigma1.*randn(nt,m);

% Plot error for several gages
figure
lw = 2;
ind = 1;
for gg=[1,5,10]
    subplot(1,3,ind)
    plot(tv, true_discharge(:,gg), 'linewidth', lw)
    hold on
    plot(tv, gage_w_error(:,gg), 'red.', 'markersize', 20)
    xlabel('Time')
    legend('Actual discharge','Measured discharge')
    ylabel('Discharge (mmd/day)')
    ind = ind + 1;
end
% ensemble ISR: multiplicative LN error with mean, stddev

err = true_discharge(:) - gage_w_error(:);
mean(err)
std(err)

figure
histogram(err)

%% Plot discharge "measurements" (Y20)

ms = 30;
fs = 16;
figure
plot(tv, true_discharge(:,1), 'linewidth', lw)
hold on 
plot(tv, gage_w_error(:,1), '.', 'markersize', ms)
legend('Truth','Measurements')
xlabel('time')
ylabel('Q (mm/day)')
grid on
title('Discharge at outlet')
set(gca, 'fontsize', fs)

% true_discharge_w_swot_sampling

%% Do ISR (Y20)

% 173 seconds per window to calculate Cx -> about 16 hours for 334 windows
% Approx. 180 seconds per window to calculate K using pinv or \ (use pinv)
% It ends up taking about 6 minutes per window to run Y20 ISR for Ohio

% runoff_init = ones(nt,n);

% s = 2*(k+1);
s = k+1;
optionsfile = './ohio_data/options_ohio.txt';
rho_thres = exp(-2);
L = 40; % km
T = 5; % days
alpha1 = 1;
R = (0.15^2);

tic
[post_runoff_Y20, Klast] = ISR_Y20(tmpa_runoff_prior, HH, gage, s, basin, optionsfile, L, T, rho_thres);
toc

% Save results
save('./ohio_data/ISR_results_Y20_swot.mat', 'post_runoff_Y20', '-v7.3')

%% Plot overview of results

truth = struct();
truth.total_runoff = nldas_runoff_true';
truth.true_runoff = nldas_runoff_true';
basin.true_runoff = truth.true_runoff;

gi = (k+1):nt-(k+1);
plt_ISR_results_overview(basin, tmpa_runoff_prior', post_runoff_Y20', truth, tv, gi)
