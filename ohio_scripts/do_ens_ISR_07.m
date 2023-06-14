%% Domain ISR with log transform and ensemble Kalman update
%
% 6/7/2023 JRS

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
t=73:77;
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

%% Assemble runoff prior (Ensemble)

% Load in the cosmos errors and use them to corrupt the true runoff

tic
M = 200; % ensemble size
A = load('./ohio_data/cosmos_error_gen/m1s1L40T5/errs/simulated_runoff_errors_200_m1s1L40T5.mat');

mean(A.alldata(:))
std(A.alldata(:))

A.alldata = A.alldata(:,:,1:nt,:);

aa=A.alldata(1,1,:,1);
figure,plot(aa(:))

figure
histogram(A.alldata(:))
xlim([0,6])
xlabel('Runoff error (mm/day)')
ylabel('Count')

runoff_errors = reshape(A.alldata, 72*104, nt, M);
basinmask = basin.mask; % no nan values
basinmask(isnan(basin.mask)) = 0;
basinmask = logical(basinmask);
runoff_errors = runoff_errors(basinmask, :, :);

% Corrupt TMPA runoff with random errors to generate prior
prior_runoff_ens = zeros(n, nt, M);
for mm=1:M
    prior_runoff_ens(:,:,mm) = tmpa_runoff_prior'.*runoff_errors(:,:,mm);
end
disp('Generated prior ensemble')
toc

% Plot prior maps

figure
t=73:77;
for i=1:5
    prior_runoff_map_ENS = make_map(basin, squeeze(mean(prior_runoff_ens(:,t(i),:),3)));
    subplot(1,5,i)
    plotraster(basin.lonv, basin.latv, prior_runoff_map_ENS, ['Prior ensemble mean (day ' num2str(t(i)) ')'])
    caxis([0,6])
end
colormap cool

%% Generate discharge "measurements" (Ensemble)

% ENKF() assumes homoscedastic, uncorrelated lognormal errors here
% we will use unbiased relative error of 15% (COV=0.15, mu = 1)
mean1 = 1;
Cv = (0.15)^2; % variance
mu1 = log((mean1^2)/sqrt(Cv+mean1^2));
sigma1 = sqrt(log(Cv/(mean1^2)+1));
error_corrupted_discharge_meas = true_discharge_w_swot_sampling.*lognrnd(mu1, sigma1, nt, m);
% error_corrupted_discharge_meas = true_discharge.*lognrnd(mu1, sigma1, nt, m);

lrv = lognrnd(mu1, sigma1, 1e4, 1);
rv = log(lrv);

% figure
% plot(tv, true_discharge)
% hold on
% plot(tv, error_corrupted_discharge_meas)
% legend('True discharge','Synthetic measurements')
% xlabel('Time')
% ylabel('Q (mm/day)')

% Plot error for several gages
figure
lw = 2;
ind = 1;
for gg=[1,5,10]
    subplot(1,3,ind)
    plot(tv, true_discharge(:,gg), 'linewidth', lw)
    hold on
    plot(tv, error_corrupted_discharge_meas(:,gg), 'red.', 'markersize', 20)
    xlabel('Time')
    legend('Actual discharge','Measured discharge')
    ylabel('Discharge (mmd/day)')
    ind = ind + 1;
end

%% Do ISR (Ensemble)

% For some reason I need to load this instance of "basin"...
% Probably due to flipud tbh
C = load('./ohio_data/swot_like_measurements_100m_no_error_revised.mat');

% Trim down to size for test purposes
% tmax = 50;
% Mmax = 10; 
% prior_runoff_ens = prior_runoff_ens(:,1:tmax,1:Mmax);
% error_corrupted_discharge_meas = error_corrupted_discharge_meas(1:tmax,:);
% % A.basin
% tv = tv(1:tmax);

% Could the issue be with make_map/flipud? Could that be why we are getting
% NaN/Inf?
opt = struct();
runoff_prior = ones(n, nt, M);
Cv = (0.15)^2;
s = k+1;
C.basin.flow_vel = flow_vel;
C.basin.timestep = 3600*24;

tic
% Takes about 8.5 hours
[post_runoff, w1, w2, maxv, maxv_sub] = ISR_domain_wrapper(prior_runoff_ens, ...
    HH, error_corrupted_discharge_meas, s, C.basin, Cv, opt, tv);
toc

% Save results
save('./ohio_data/ISR_results_domain_m240_swot.mat', 'post_runoff', ...
    'error_corrupted_discharge_meas','s','Cv','opt','tv','-v7.3')

% post_runoff_orig = post_runoff;
% post_runoff(post_runoff>1e5) = 0;
% post_runoff(isnan(post_runoff)) = 0;

% post_runoff_orig = post_runoff;
% imag_ind = find(abs(imag(post_runoff))>0);
% post_runoff(imag_ind) = NaN;

%%

ENS = struct(); % structure to store results
ENS.mean_prior_runoff = mean(prior_runoff_ens,3);       
ENS.sd_prior_runoff = std(prior_runoff_ens, [], 3);
ENS.mean_posterior_runoff = mean(post_runoff, 3);
ENS.sd_posterior_runoff = std(post_runoff, [], 3);

truth = struct();
truth.total_runoff = nldas_runoff_true';
truth.true_runoff = nldas_runoff_true';
basin.true_runoff = truth.true_runoff;
gi = (k+1):nt-(k+1);
% plot_ensemble_runoff_maps(basin, ENS.mean_prior_runoff, ENS.mean_posterior_runoff, truth, tv, gi)
plt_ISR_results_overview(basin, ENS.mean_prior_runoff, ENS.mean_posterior_runoff, truth, tv, gi)

