% Hetch Hetchy ISR for paper
%
% 6/14/2023
% Results don't look super good for Ensemble ISR...

clear, clc, close all
addpath(genpath('./src/'))

load('./hh_data/isr_setup_5.mat')

% save('./Data/isr_setup1.mat', 'basin', 'truth', 'flow_vel', 'k', 'HH', 'fdir', 'gage')

% basin.distmat = calc_distance_matrix(basin.lon, basin.lat);
% basin.true_runoff = truth.total_runoff';

ms = 20;
figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'HH2')
hold on
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)

%% Cut to a particular time period

timesteps = 5000:5500;
% timesteps = 5000:7000;
true_discharge = truth.discharge_mm(timesteps,:); % water year 2006
truth.total_runoff = truth.total_runoff(:,timesteps);
tv = truth.t(timesteps);

%% Plot discharge at each gauge

m = 5; % number of gages
figure
for i=1:m
    plot(tv, true_discharge(:,i))
    hold on
end
legend('g1','g2','g3','g4','g5')

%% Make "SWOT sampling"

true_discharge_w_swot_sampling = true_discharge;

%% Corrupt true discharge with error

gage_all = true_discharge;
gage = true_discharge; % should corrupt with error
% gage = true_discharge_w_swot_sampling; % should corrupt with error

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
for gg=[1,3,5]
    subplot(1,3,ind)
    plot(tv, gage_all(:,gg), 'linewidth', lw)
    hold on
    plot(tv, gage_w_error(:,gg), 'red.', 'markersize', 20)
    xlabel('Time')
    legend('Actual discharge','Measured discharge')
    ylabel('Discharge (mmd/day)')
    ind = ind + 1;
end
% ensemble ISR: multiplicative LN error with mean, stddev

err = true_discharge(:) - gage_w_error(:);
nanmean(err)
nanstd(err)

figure
histogram(err)

%% Run ISR

[n,nt] = size(truth.total_runoff);
s = 2*(k+1);
% s = k+18;
total_mean = mean(truth.total_runoff(:)); % mm/hour

alpha = 1;
% CV = (10^2)*eye(m*(s+1));
CV = 0.15^2; % meas error covariance
% CV = 0.15^2; % meas error covariance

runoff_init = total_mean*ones(n,nt);
[post_runoff,~,w1,w2] = ISR_PW13(runoff_init', HH, gage_w_error, s, 'proportional', alpha^2, CV);
post_runoff = post_runoff';

% confine analysis to a single assimilation window
i1=w1(1);
i2=w2(1);

lw = 2;
fs = 14;

%% Basin average runoff (PW13 method)

saveflag = 0;

figure
h1 = plot(tv, mean(runoff_init,1), 'b', 'linewidth', lw);
hold on
% plot(tv, mean(runoff_init,1) + mean(sd_prior_runoff,1), '--b', 'linewidth', lw)
% plot(tv, mean(runoff_init,1) - mean(sd_prior_runoff,1), '--b', 'linewidth', lw)
h2 = plot(tv, mean(post_runoff,1), 'r', 'linewidth', lw);
% plot(tv, mean(mean_posterior_runoff,1) + mean(sd_posterior_runoff,1), 'r--', 'linewidth', lw)
% plot(tv, mean(mean_posterior_runoff,1) - mean(sd_posterior_runoff,1), 'r--', 'linewidth', lw)
h3 = plot(tv, mean(truth.total_runoff,1), 'linewidth', lw, 'color', 'black');
grid on
% xlim([tv(i1), tv(i2)])
xlabel('Time')
ylabel('Runoff (mm/day)')
title('Basin-mean runoff for a single window')
set(gca, 'fontsize', fs)
if saveflag
    saveas(gcf, fullfile(wkdir, ['basin_avg_runoff_outputs.png']))
end
% draw vertical lines to indicate full vs. partial update
xline(tv(i2-k), 'black')
% the last k steps aren't fully updated in a window
% however, these results have been blended with the full update from the
% next time window, so you can't see the effect of the partial update
legend([h1,h2,h3], {'Prior Mean','Posterior Mean','Truth'})

%% Generate prior ensemble of runoff for ISR (ens)

rng(123, 'twister')

% Generate prior ensemble
mean1 = total_mean;
stddev = alpha*mean1;
mu1 = log((mean1^2)/sqrt(stddev^2+mean1^2));
sigma1 = sqrt(log(stddev^2/(mean1^2)+1));
% for mm=1:1000
%     w = lognrnd(mu1, sigma1, 1e5, 1);
%     mean(w);
%     stddevs(mm) = std(w);
%     mm;
% end
% figure,histogram(stddevs)
% mean(stddevs)
% when you generate lognormal random variables with high standard
% deviations, the sample statistics get thrown off by outliers

M = 500;
w = lognrnd(mu1, sigma1, n, nt, M);
mean(w(:))
std(w(:))
% Postulate the discharge measurement error covariance matrix
prior_runoff_ens = w;

%% Do ensemble ISR

rng(704753262, 'twister')

basin.mask(isnan(basin.mask)) = 0;
basin.mask = logical(basin.mask);
tic
Cv = 0.15^2;
options.loc = 0;
options.locL = 0.1; % cells
options.locT = 0.1; % hr
basin.distmat = calc_distance_matrix(basin.grid_lon(basin.mask), basin.grid_lat(basin.mask));
basin.true_runoff = truth.total_runoff;
basin.t = tv;
basin.true_discharge = true_discharge;
[post_runoff, w1, w2] = ISR_EnKF(prior_runoff_ens, HH, true_discharge, s, basin, Cv, options);
toc

mean_prior_runoff = mean(prior_runoff_ens,3);       
sd_prior_runoff = std(prior_runoff_ens, [], 3);
mean_posterior_runoff = mean(post_runoff, 3);
sd_posterior_runoff = std(post_runoff, [], 3);

%% Also do ISR with domain localization

basin.flow_vel = 2;
basin.timestep = 3600;
[post_runoff, w1, w2] = ISR_domain_wrapper(prior_runoff_ens, HH, true_discharge, s, basin, Cv, options, tv);

%% Basin average runoff (ensemble)

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
xline(tv(i2-k), 'black')
% the last k steps aren't fully updated in a window
% however, these results have been blended with the full update from the
% next time window, so you can't see the effect of the partial update
legend([h1,h2,h3], {'Prior Mean','Posterior Mean','Truth'})
ylim([-10,10])

%% Look at predicted measurements at each gage

mean_prior_discharge = state_model_dumb(mean_prior_runoff', HH);
mean_posterior_discharge = state_model_dumb(mean_posterior_runoff', HH);

m = 5;
prior_discharge = zeros(nt,m,M);
posterior_discharge = zeros(nt,m,M);
for mm=1:M
    prior_discharge(:,:,mm) = state_model_dumb(prior_runoff_ens(:,:,mm)', HH);
    posterior_discharge(:,:,mm) = state_model_dumb(post_runoff(:,:,mm)', HH);
end

prior_discharge_sd = std(prior_discharge, [], 3);
posterior_discharge_sd = std(posterior_discharge, [], 3);

%% Gage predictions

figure

for g=1:5
    
subplot(2,3,g)    
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
%     h4 = plot(tv, error_corrupted_discharge_meas(:,g), '.green', 'markersize', ms);
    
%     xlim([tv(i1), tv(i2)])
%     xlim([datetime(2006,5,13), datetime(2006,5,14)])
%     xlim([datetime(2006,5,13), datetime(2006,5,18)])
    title(['Gage ' num2str(g)])
    xlabel('Date')
    ylabel('Discharge (mmd)')
    if g==11
        legend([h1,h2,h3,h4], {'Prior mean','Posterior mean','Truth', 'Measurements'})
    end
    if saveflag 
        saveas(gcf, fullfile(wkdir, ['discharge_at_gage_' num2str(g) '.png']))
    end
    
end

