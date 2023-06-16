% Hetch Hetchy ISR for paper
%
% 12/19/2022

clear, clc, close all
addpath(genpath('./src/'))

% create inputs using basin_prep_HH2.m
% load('/hdd/ISR/02_Basins/HH2/Data/isr_data.mat')
load('./hh_data/isr_setup_5.mat')
% load('./Data/unsorted/isr_setup4.mat')
% load('./Data/basin_setups/isr_setup_usds_gage_2.mat')
% load('./Data/basin_setups/isr_setup_us_gage_1.mat')

% load('/hdd/ISR/02_Basins/HH2/Data/unsorted/isr_setup1.mat')

% basin.lonv = lonv;
% basin.latv = latv;
% save('./Data/isr_setup1.mat', 'basin', 'truth', 'flow_vel', 'k', 'HH', 'fdir', 'gage')

% basin.distmat = calc_distance_matrix(basin.lon, basin.lat);
% basin.true_runoff = truth.total_runoff';

ms = 20;
% figure
% plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'HH2')
% hold on
% plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)

basin.mask(isnan(basin.mask)) = 0;
basin.mask = logical(basin.mask);

%% Cut to a particular time period

i1 = find(truth.t == datetime(2006,5,1));
i2 = find(truth.t == datetime(2006,5,31));
timesteps = i1:i2; % originally 5000:7000

% looking at April and May 2007

% timesteps = 5000:7000;
true_discharge = truth.discharge_mm(timesteps,:); % water year 2006
truth.total_runoff = truth.total_runoff(:,timesteps);
tv = truth.t(timesteps);

% cropping to the time period from Oct. 1, 2006 to September 30, 2007

%% Plot discharge at each gauge

fs = 18;
m = size(true_discharge,2); % number of gages
figure
for i=1:m
    plot(tv, true_discharge(:,i))
    hold on
end
xlabel('Time')
ylabel('Q (mm/day)')
title('True discharge at Hetch Hetchy gages')
legend('g1','g2','g3','g4','g5')
set(gca, 'fontsize',fs)

%% EnKF implementation

s = 2*k+1;
mn = 1;
sd = 0.15;
error_corrupted_discharge_meas = generate_discharge_measurements(true_discharge, 'norm', mn, sd, 1,0);
error_corrupted_discharge_meas_swot = generate_discharge_measurements(true_discharge, 'norm', mn, sd, 1,1);

%% Use mvlognrnd to generate uncorrelated errors

addpath('/home/jschap/Documents/MATLAB/MvLogNRand')
rng(123, 'twister')

% case 2
% mean1 = 1;
% stddev = 0;

% case 4
mean1 = 2;
stddev = 2;

mu1 = log((mean1^2)/sqrt(stddev^2+mean1^2));
sigma1 = sqrt(log(stddev^2/(mean1^2)+1));

M = 500;
% w = lognrnd(mu1, sigma1, n, nt, M);

CorrMat = eye(n);
w = zeros(n,nt,M);
for mm=1:M
    w1 = MvLogNRand(mu1*ones(n,1), sigma1*ones(n,1), nt, CorrMat);
    w(:,:,mm) = w1';
end

%% Generate runoff prior (uncorr)

prior_runoff_ens = zeros(n, nt, M);
for m=1:M
    prior_runoff_ens(:,:,m) = truth.total_runoff.*w(:,:,m);
end

%% Generate runoff prior (perfectly correlated, biased)

prior_runoff_ens = zeros(n, nt, M);
w = 2*ones(n,nt,M);
for m=1:M
    prior_runoff_ens(:,:,m) = truth.total_runoff.*w(:,:,m);
end

%% Load ensemble of runoff errors

% wkdir = '/hdd/ISR/02_Basins/HH2/Data/prior_err_L1_S24_mu2_sigma_2';
wkdir = '/hdd/ISR/02_Basins/HH2/Data/exp1-4/exp3';
% wkdir = '/hdd/ISR/02_Basins/HH2/Data/runoff_priors/L538_T117_mu2_sigma2';
% wkdir = '/hdd/ISR/02_Basins/HH2/sensitivity_to_prior/varyingL/L10000/';
A = load(fullfile(wkdir, 'simulated_runoff_errors_500.mat'));

M = 500; % ensemble size

mean(A.alldata(:))
std(A.alldata(:))

A.alldata = A.alldata(:,:,1:nt,:);

aa=A.alldata(1,1,:,1);
figure,plot(aa(:))

autocorr(aa(:), 'numlags', 200)

figure
histogram(A.alldata(:))
xlim([0,6])
xlabel('Runoff error (mm/day)')
ylabel('Count')

runoff_errors = reshape(A.alldata, 6*7, nt, M);
runoff_errors = runoff_errors(basin.mask, :, :);

prior_runoff_ens = zeros(n, nt, M);
for m=1:M
    prior_runoff_ens(:,:,m) = truth.total_runoff.*runoff_errors(:,:,m);
end

%% Plot prior runoff ensemble

basin.lon=basin.lonv; basin.lat=basin.latv;
% basin.mask = flipud(basin.mask);

lw = 2;
tmp = squeeze(mean(prior_runoff_ens,1));

figure
h1 = plot(tv, tmp, 'cyan');
hold on
h2 = plot(tv, mean(tmp,2), '--*b', 'linewidth', lw);
h3 = plot(tv, mean(truth.total_runoff,1), 'k', 'linewidth', lw);
xlabel('Time')
ylabel('Runoff (mm/hour)')
legend([h1(1);h2;h3], {'Ensemble','Ens mean','Truth'})
title('Basin average runoff')
xlim([datetime(2006,5,13), datetime(2006,5,16)])
set(gca, 'fontsize', fs)
% saveas(gcf, fullfile(wkdir, ['basin_avg_runoff_ensemble_inputs.png']))

% Plot prior runoff ensemble for a particular grid cell (cell 2/outlet)

% kk = 2; % grid cell
% figure(8)
% subplot(2,2,4)
% h1 = plot(tv, squeeze(prior_runoff_ens(kk,:,:)), 'cyan');
% hold on
% h2 = plot(tv, mean(prior_runoff_ens(kk,:,:),3), '--*b', 'linewidth', lw);
% h3 = plot(tv, truth.total_runoff(kk,:), 'k', 'linewidth', lw);
% xlim([datetime(2006,5,13), datetime(2006,5,18)])
% xlabel('Time')
% ylabel('Runoff (mm/day)')
% % legend([h1(1);h2;h3], {'Ensemble','Ens mean','Truth'})
% title('Scenario 4')
% set(gca, 'fontsize', fs)

% Plot map of prior ensemble
figure
aa = squeeze(prior_runoff_ens(:,20,:));
% aa = squeeze(mean(prior_runoff_ens,2));
for kk=1:9
    subplot(3,3,kk)
    aaa = make_map(basin, aa(:,kk));    
    plotraster(basin.lon, basin.lat, aaa, ['Ensemble member num2str(kk)'])
end

%% Estimate localization parameters from the runoff prior

% Could theoretically be done, I think. But we know them a priori, so we
% don't need to estimate them

% These are Cyy parameters. We assume they also apply for Cyz and Czz.

% We could try using C. David's (2019) method for propagating runoff error to
% discharge error, but that is out of scope

%% Run EnKF-ISR (RUN)

rng(704753262, 'twister')

% Postulate the discharge measurement error covariance matrix
% basin.mask = flipud(basin.mask);

s=k+1;
s=30;
tic
options.loc = 0;
options.locL = 5; % cells
options.locT = 72; % hr
options.rho_yy=0;
basin.distmat = calc_distance_matrix(basin.grid_lon(basin.mask), basin.grid_lat(basin.mask));
basin.true_runoff = truth.total_runoff;
basin.t = tv;
basin.true_discharge = true_discharge;

% with discharge error and SWOT sampling
% [post_runoff, w1, w2] = ISR_EnKF(prior_runoff_ens, HH, error_corrupted_discharge_meas_SWOT, s, basin, Cv, options);
% [post_runoff, w1, w2, maxv, maxv_sub] = ISR_EnKF_domain_wrapper_d2(prior_runoff_ens, HH, error_corrupted_discharge_meas_SWOT, s, basin, Cv, options, tv);

% with discharge error
[post_runoff, w1, w2] = ISR_EnKF(prior_runoff_ens, HH, error_corrupted_discharge_meas, s, basin, Cv, options);
% [post_runoff, w1, w2, maxv, maxv_sub] = ISR_EnKF_domain_wrapper_d2(prior_runoff_ens, HH, error_corrupted_discharge_meas, s, basin, Cv, options, tv);

% no discharge error
% [post_runoff, w1, w2] = ISR_EnKF(prior_runoff_ens, HH, true_discharge, s, basin, 0.001^2, options);
% [post_runoff, w1, w2] = ISR_EnKF_domain_wrapper_d2(prior_runoff_ens, HH, true_discharge, s, basin, 0.001^2, options, tv);

toc

mean_prior_runoff = mean(prior_runoff_ens,3);       
sd_prior_runoff = std(prior_runoff_ens, [], 3);
mean_posterior_runoff = mean(post_runoff, 3);
sd_posterior_runoff = std(post_runoff, [], 3);

saveflag = 0;

fs = 18;
lw = 2;

% confine analysis to a single assimilation window
i1=w1(1);
i2=w2(1);
% xlim([tv(i1), tv(i2)])
% don't look at the spin-up/cooldown part of the window

% average bias of the prior
mean(prior_runoff_ens(:))/mean(truth.total_runoff(:))

% average bias of the posterior
mean(post_runoff(:))/mean(truth.total_runoff(:))

% seems to be a -4.24% bias in the posterior runoff. Why? Shouldn't be
% any... the method seems to introduce it. Can we get rid of it?

%% Y20 method

L = 5;
T = 6;
rho_thres = exp(-2);
runoff_init = truth.total_runoff;
optionsfile = './hh_data/options_hh2.txt';
[runoff, K] = ISR_Y20(runoff_init', HH, error_corrupted_discharge_meas, ...
    s, basin, optionsfile, L, T, rho_thres); % Y20 method

%% Maps

gi = 1:nt;
basin.mask = flipud(basin.mask);
plt_ISR_results_overview(basin, mean_prior_runoff, mean_posterior_runoff, truth, tv, gi)
% plot_ensemble_runoff_maps(basin, mean_prior_runoff, mean_posterior_runoff, truth, tv, 1:nt, 1)
% plot_ensemble_runoff_maps(basin, B.mean_prior_runoff, B.mean_posterior_runoff, truth, tv)
% plot_ensemble_runoff_maps(basin, runoff_init, runoff', truth, tv, 1:nt, 1)

%% Figure for paper

gi=1:nt;
ms=15;

figure
subplot(1,2,1)
plot_gofmaps(basin, B.mean_posterior_runoff', truth.total_runoff, gi);
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = 'NSE';
set(colorTitleHandle ,'String',titleString);
title('Without localization')
hold on
plot(basin.gage_lon, basin.gage_lat,'bo','markersize',ms,'linewidth',lw)

subplot(1,2,2)
plot_gofmaps(basin, mean_posterior_runoff', truth.total_runoff, gi);
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = 'NSE';
set(colorTitleHandle ,'String',titleString);
title('With localization')
hold on
plot(basin.gage_lon, basin.gage_lat,'bo','markersize',ms,'linewidth',lw)

%% Save results

save('./effect_of_domain_loc/usds_gage_2_results_L1_T24.mat')
B = load('./effect_of_domain_loc/usds_gage_2_results_L1_T24.mat');

%% Look at EnKF ISR outputs (RUN)

figure
h1 = plot(truth.t(timesteps), mean(mean_prior_runoff,1), 'b', 'linewidth', lw);
hold on
plot(truth.t(timesteps), mean(mean_prior_runoff,1) + mean(sd_prior_runoff,1), '--b', 'linewidth', lw)
plot(truth.t(timesteps), mean(mean_prior_runoff,1) - mean(sd_prior_runoff,1), '--b', 'linewidth', lw)
h2 = plot(truth.t(timesteps), mean(mean_posterior_runoff,1), 'r', 'linewidth', lw);
plot(truth.t(timesteps), mean(mean_posterior_runoff,1) + mean(sd_posterior_runoff,1), 'r--', 'linewidth', lw)
plot(truth.t(timesteps), mean(mean_posterior_runoff,1) - mean(sd_posterior_runoff,1), 'r--', 'linewidth', lw)
h3 = plot(truth.t(timesteps), mean(truth.total_runoff,1), 'linewidth', lw, 'color', 'black');
grid on
% ylim([-1,20])
% xlim([tv(i1), tv(i2)])
% xlim([datetime(2006,5,13), datetime(2006,5,14)])
% xlim([datetime(2006,5,13), datetime(2006,5,18)])
xlabel('Time')
ylabel('Runoff (mm/day)')
title('Basin-mean runoff for a single window')
set(gca, 'fontsize', fs)
if saveflag
    saveas(gcf, fullfile(wkdir, ['basin_avg_runoff_outputs.png']))
end
% draw vertical lines to indicate full vs. partial update
% xline(tv(i2-k), 'black')
% the last k steps aren't fully updated in a window
% however, these results have been blended with the full update from the
% next time window, so you can't see the effect of the partial update
legend([h1,h2,h3], {'Prior Mean','Posterior Mean','Truth'})

%% Look at predicted measurements at each gage

mean_prior_discharge = state_model_dumb(mean_prior_runoff', HH);
mean_posterior_discharge = state_model_dumb(mean_posterior_runoff', HH);

% m = 5;
prior_discharge = zeros(nt,m,M);
posterior_discharge = zeros(nt,m,M);
for mm=1:M
    prior_discharge(:,:,mm) = state_model_dumb(prior_runoff_ens(:,:,mm)', HH);
    posterior_discharge(:,:,mm) = state_model_dumb(post_runoff(:,:,mm)', HH);
end

prior_discharge_sd = std(prior_discharge, [], 3);
posterior_discharge_sd = std(posterior_discharge, [], 3);

%%
% figure

for g=1:m
    
    figure
% subplot(2,3,g)    
%     figure(8)
%     subplot(2,1,2)
    h1 = plot(tv, mean_prior_discharge(:,g), 'b', 'linewidth', lw);
    hold on
    plot(tv, mean_prior_discharge(:,g) - prior_discharge_sd(:,g), 'b--', 'linewidth', lw);
    plot(tv, mean_prior_discharge(:,g) + prior_discharge_sd(:,g), 'b--', 'linewidth', lw);
    
    h2 = plot(tv, mean_posterior_discharge(:,g), 'r-+', 'linewidth', lw, 'markersize', 10);
    plot(tv, mean_posterior_discharge(:,g) - posterior_discharge_sd(:,g), 'r--', 'linewidth', lw);
    plot(tv, mean_posterior_discharge(:,g) + posterior_discharge_sd(:,g), 'r--', 'linewidth', lw);
    
    h3 = plot(tv, true_discharge(:,g), 'k', 'linewidth', lw);
    h4 = plot(tv, error_corrupted_discharge_meas(:,g), '.green', 'markersize', ms);
    
%     xlim([tv(i1), tv(i2)])
%     xlim([datetime(2006,5,13), datetime(2006,5,14)])
%     xlim([datetime(2006,5,6), datetime(2006,5,9)])
    title(['Gage ' num2str(g)])
    xlabel('Date')
    ylabel('Discharge (mmd)')
    if g==1
        legend([h1,h2,h3,h4], {'Prior mean','Posterior mean','Truth', 'Measurements'})
    end
    if saveflag 
        saveas(gcf, fullfile(wkdir, ['discharge_at_gage_' num2str(g) '.png']))
    end
    
end

%% Figure for Steve (one time window/deep dive into update/discharge)

g = 1;
figure
ttv = tv(i1:i2);
    
h1 = plot(mean_prior_discharge(i1:i2,g), 'b', 'linewidth', lw);
hold on
plot(mean_prior_discharge(i1:i2,g) - prior_discharge_sd(i1:i2,g), 'b--', 'linewidth', lw);
plot(mean_prior_discharge(i1:i2,g) + prior_discharge_sd(i1:i2,g), 'b--', 'linewidth', lw);

h2 = plot(mean_posterior_discharge(i1:i2,g), 'r', 'linewidth', lw);
plot(mean_posterior_discharge(i1:i2,g) - posterior_discharge_sd(i1:i2,g), 'r--', 'linewidth', lw);
plot(mean_posterior_discharge(i1:i2,g) + posterior_discharge_sd(i1:i2,g), 'r--', 'linewidth', lw);

h3 = plot(true_discharge(i1:i2,g), 'k', 'linewidth', lw);
h4 = plot(error_corrupted_discharge_meas(i1:i2,g), '.g', 'markersize', ms);
title(['Gage ' num2str(g)])
ylabel('Discharge (mmd)')
legend([h1,h2,h3,h4], {'Prior mean','Posterior mean','Truth', 'Measurements'})
grid on
set(gca, 'fontsize', fs)
if saveflag
    saveas(gcf, fullfile(wkdir, ['single_update_window.png']))
end

%% Look at map of prior and posterior runoff NSE (RUN)

gi = 7*24:nt-k-1;

figure
subplot(2,1,1)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, mean_prior_runoff', truth.total_runoff, gi);
title('Prior runoff NSE')
subplot(2,1,2)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, mean_posterior_runoff', truth.total_runoff, gi);
title('Posterior runoff NSE')
if saveflag
    saveas(gcf, fullfile(wkdir, ['gofmaps.png']))
end

%% Look at ensemble of runoff at individual grid cells

kk1 = [2,7,12,17]; % cell index
lw2 = 0.5;

figure(1)
for mm=1:4
    
    subplot(2,4,mm)
    
    kk=kk1(mm); % indexing
    a1 = plot(tv, squeeze(prior_runoff_ens(kk, :, :)), 'cyan', 'linewidth', lw2);
    hold on
    a2 = plot(tv, squeeze(post_runoff(kk, :, :)), 'magenta', 'linewidth', lw2);
    a3 = plot(tv, mean_prior_runoff(kk,:), 'blue', 'linewidth', lw);
    a4 = plot(tv, mean_posterior_runoff(kk,:), 'green', 'linewidth', lw);
    a5 = plot(tv, truth.total_runoff(kk,:), 'black', 'linewidth', lw);
%     legend([a1(1),a2(1),a3,a4, a5], {'Prior ensemble','Posterior ensemble','Prior mean','Posterior mean', 'Truth'}, ...
%         'location', 'best')
    grid on
    title(['Grid cell ' num2str(kk) ' runoff (mmd)'])
    % set(gca, 'fontsize', fs)
    % xlim([ttv(1), ttv(end)])
%     ylim([0,8])
    xlim([datetime(2006,5,13), datetime(2006,5,14)])
    if saveflag
        saveas(gcf, fullfile(wkdir, ['ISR_results_for_grid_cell_' num2str(kk) '.png']))
    end

end

%% Figure for paper (look at individual grid cells)

kk1=1:20;
lw2 = 0.5;

% Find cells with largest posterior error
sdposterior = zeros(n,1);
for mm=1:n
    xposterior = squeeze(post_runoff(mm, :, :));
    sdposterior(mm) = mean(std(xposterior,[],2));
end
[~,minind] = min(sdposterior);
[~,maxind] = max(sdposterior);

kk1 = [minind, maxind]; % cell index

kk1 = [2,11,17];

figure(1)
for mm=1:length(kk1)
    
    subplot(1,length(kk1),mm)
    
    kk=kk1(mm); % indexing
        
    xprior = squeeze(prior_runoff_ens(kk, :, :));
    xposterior = squeeze(post_runoff(kk, :, :));
    sdprior = std(xprior,[],2);
    sdposterior = std(xposterior,[],2);

    a1 = plot(tv, mean_prior_runoff(kk,:)' - sdprior, 'blue--', 'linewidth', lw);
    hold on
    a2 = plot(tv, mean_prior_runoff(kk,:)' + sdprior, 'blue--', 'linewidth', lw);
    a3 = plot(tv, mean_posterior_runoff(kk,:)' - sdposterior, 'red--', 'linewidth', lw);
    a4 = plot(tv, mean_posterior_runoff(kk,:)' + sdposterior, 'red--', 'linewidth', lw);
    a7 = plot(tv, truth.total_runoff(kk,:), 'black-+', 'linewidth', lw);
    a5 = plot(tv, mean_prior_runoff(kk,:), 'blue', 'linewidth', lw);
    a6 = plot(tv, mean_posterior_runoff(kk,:), 'red', 'linewidth', lw);
    
    if mm==1
    legend([a5,a6,a7], {'Prior mean','Posterior mean','Truth'}, ...
        'location', 'best')
    end
    
    grid on
    title(['Grid cell ' num2str(kk)])
    % set(gca, 'fontsize', fs)
    xlim([tv(90), tv(210)])
%     ylim([0,8])
    xlabel('Time')
    ylabel('Runoff (mm/hour)')
    set(gca, 'fontsize', fs)
%     xlim([datetime(2006,5,13), datetime(2006,5,14)])
%     xlim([tv(i1), tv(i2)])

end

% find(tv==datetime(2006,5,1))
% find(tv==datetime(2006,5,6))

%% Map of runoff (snapshot in time)

[nr,nc] = size(basin.mask);
maps = struct();
maps.true_runoff = zeros(nr,nc,nt);
maps.mean_prior_runoff = zeros(nr,nc,nt);
maps.mean_post_runoff = zeros(nr,nc,nt);
for t=1:nt
    maps.true_runoff(:,:,t) = make_map(basin, truth.total_runoff(:,t));
    maps.mean_prior_runoff(:,:,t) = make_map(basin, mean_prior_runoff(:,t));
    maps.mean_post_runoff(:,:,t) = make_map(basin, mean_posterior_runoff(:,t));
end
    
cmax = 6;
i1 = 64;
figure
subplot(1,3,1)
plotraster('lon','lat',maps.true_runoff(:,:,i1), ['Truth ' datestr(tv(i1))])
caxis([0,cmax])
set(gca, 'fontsize', 12)
subplot(1,3,2)
plotraster('lon','lat',maps.mean_prior_runoff(:,:,i1), 'Prior mean (mm/hr)')
caxis([0,cmax])
subplot(1,3,3)
plotraster('lon','lat',maps.mean_post_runoff(:,:,i1), 'Posterior mean')
caxis([0,cmax])

%% Map of upstream cells (downstream cells?)

basin.ds_gages = flipud(sum(basin.mask_gage,3));
figure
plotraster('lon','lat',basin.ds_gages, 'downstream gages')

% Flow accumulation (number of cells) is a predictor of how good the update
% is. This is because there are more cells to gauges, more unknowns to
% knowns.

%% Map of basin with cell numbers and gages

n = size(HH,2);
[nr, nc] = size(basin.mask);
basin.cellnum = double(basin.mask);
basin.cellnum(basin.mask) = (1:n)';

x = repmat(1:nc,nr,1); % generate x-coordinates
y = repmat((1:nr)',1,nc); % generate y-coordinates
% Generate Labels
t = num2cell(basin.cellnum); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% Draw Image and Label Pixels
figure
imagesc(basin.mask)
text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
title('Grid cell numbers')

%% Reduce gages and preconditioning for domain ISR (not needed)

% This function checks the gauge network to see if there will be negative
% discharges and suggests which gauges to remove to enable domain ISR
rm_ind = domain_isr_precond(basin, error_corrupted_discharge_meas);

isr_setup_in = './Data/unsorted/isr_setup5.mat';

isr_setup_out = './Data/unsorted/isr_setup4.mat';
reduce_gages(isr_setup_in, isr_setup_out, rm_ind); % in this case, we remove gauge 2, since it causes problems for domain ISR (it is too close to gage 1)

% save('./Data/isr_setup1.mat', 'basin', 'truth', 'flow_vel', 'k', 'HH', 'fdir', 'gage')

% Tabling this idea for now (4/21/23)
% Instead, we're going to try to avoid having those negative adjusted
% discharge values in the first place, by not adjusting the discharge
% values
%
% Instead, the idea is to fix the already estimated runoff at its current
% values, and use that in the measurement model. This will require some
% thoughtful application of the routing model