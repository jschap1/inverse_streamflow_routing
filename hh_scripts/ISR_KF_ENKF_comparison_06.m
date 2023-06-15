% Comparison between PW13 method and the new ensemble ISR method
%
% 6/14/2023
% This script can pretty much just be run without babysitting
% It generates PW13 vs. ensemble ISR figures for Hetch Hetchy and shows
% where our method is better. It is up to us to describe exactly why it is
% better, and corroborate our numerical findings with theory.

clear, clc, close all
addpath(genpath('./src/'))

load('./hh_data/isr_setup_5.mat')

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

%% Run ISR

[n,nt] = size(truth.total_runoff);
s = 2*(k+1);
% s = k+18;
total_mean = mean(truth.total_runoff(:)); % mm/hour

R = 0.15; % relative error standard deviation for meas errors (% error)
alpha1 = 1; % runoff error proportionality constant

for tt=1:nt
    Cv = eye(m);
    Cv(1:m+1:end) = (R*true_discharge(tt,:)).^2;
    v = mvnrnd(zeros(m,1), Cv);
    error_corrupted_discharge(tt,:) = true_discharge(tt,:) + v;
end

g=1; % gage number
ms=10;
figure
plot(tv, error_corrupted_discharge(:,g), 'cyan.', 'markersize', ms)
hold on
plot(tv, true_discharge(:,g), 'k')
legend('measurement','true')
xlabel('Time')
ylabel('Discharge (mm/hr)')

runoff_init = 2*truth.total_runoff;
% runoff_init = total_mean*ones(n,nt);
[post_runoff,~,w1,w2] = ISR_PW13(runoff_init', HH, error_corrupted_discharge, s, 'proportional', alpha1, R^2);
post_runoff = post_runoff';

pw13.post_runoff = post_runoff;
pw13.prior_runoff = runoff_init;

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

%% Look at posterior discharge (KF)

y_prior = state_model_dumb(runoff_init', HH); 
y_post = state_model_dumb(post_runoff', HH);

ms=10;
for g=1:5
    
	figure  
    h1 = plot(tv, y_prior(:,g), 'b', 'linewidth', lw);
    hold on   
    h2 = plot(tv, y_post(:,g), 'r-+', 'linewidth', lw, 'markersize', ms);
    h3 = plot(tv, true_discharge(:,g), 'k', 'linewidth', lw);
    h4 = plot(tv, error_corrupted_discharge(:,g), 'green.', 'markersize', ms);
%     xlim([tv(i1), tv(i2)])
    title(['Gage ' num2str(g)])
    xlabel('Date')
    ylabel('Discharge (mmd)')

end

%% Maps

gi=1:nt;

% basin.mask = flipud(basin.mask);
basin.true_runoff = truth.total_runoff;
plt_ISR_results_overview(basin, runoff_init, post_runoff, truth, tv, gi)

figure

subplot(1,2,1)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, runoff_init', truth.total_runoff, gi);
title('Prior mean (NSE)')
subplot(1,2,2)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, post_runoff', truth.total_runoff, gi);
title('Posterior mean (NSE)')

min(post_runoff(:))

i1=378;
i2=i1+24;
figure
for c = 1:n
    subplot(4,5,c)
    plot(tv, runoff_init(c,:), 'b')
    hold on
    plot(tv, post_runoff(c,:),'r')
    plot(tv, truth.total_runoff(c,:), 'k')
    if c==1
        legend('prior','posterior','truth')
    end
%     xlim([tv(i1),tv(i2)])
    title(['Runoff for cell ' num2str(c)])
end

for c=1:n
kge(c) = myKGE(truth.total_runoff(c,:)', post_runoff(c,:)');
nse(c) = myNSE(truth.total_runoff(c,:)', post_runoff(c,:)');
end

%% Maps for specific cells, different values of alpha

figure(11)
ii=1; % 1, 4, 7
i1 = find(tv==datetime(2006,5,1));
i2 = find(tv==datetime(2006,5,6));
for c=[2,11,17]
    subplot(1,3,ii)
    plot(tv(i1:i2), runoff_init(c,i1:i2), 'b','linewidth',lw)
    hold on
    plot(tv(i1:i2), post_runoff(c,i1:i2),'r','linewidth',lw)
    plot(tv(i1:i2), truth.total_runoff(c,i1:i2), 'k','linewidth',lw)    
    ii = ii+1;
    title(['Cell ' num2str(c) ', \alpha = ' num2str(alpha1)])
    xlabel('Time')
    ylabel('Runoff (mm/hr)')
    switch c
        case 2
            ylim([-1, 2.5]);
            if ii==1, legend('prior','posterior','truth'), end
        case 11
            ylim([-1, 3]);
        case 17
            ylim([0, 1.1]);
    end
    set(gca, 'fontsize', fs)
end

myKGE(truth.total_runoff(c,i1:i2)', post_runoff(c,i1:i2)')
myKGE(truth.total_runoff(c,i1:i2)', runoff_init(c,i1:i2)')

i1 = find(tv==datetime(2006,5,1));
i2 = find(tv==datetime(2006,5,6));
figure(4)
g=1; % gage at the outlet
subplot(2,1,1)
plot(tv(i1:i2), y_prior(i1:i2,g), 'b', 'linewidth', lw);
hold on   
plot(tv(i1:i2), y_post(i1:i2,g), 'r-+', 'linewidth', lw, 'markersize', ms);
plot(tv(i1:i2), true_discharge(i1:i2,g), 'k', 'linewidth', lw);
plot(tv(i1:i2), error_corrupted_discharge(i1:i2,g), 'green.', 'markersize', ms);
legend('Prior','Posterior','Truth','Meas.')
title(['\alpha = ' num2str(alpha1)])
xlabel('Time')
ylabel('Discharge(mm/hr)')
set(gca, 'fontsize', fs)

%% Generate prior ensemble of runoff for ISR (ens)

rng(123, 'twister')

% Generate prior ensemble
mean1 = 2;
stddev = 1; % this is a standard deviation
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

% Simulate uncorrelated runoff errors and create the ensemble
CorrMat = eye(n);
prior_runoff_ens = zeros(n,nt,M);
for mm=1:M
    w = MvLogNRand(mu1*ones(n,1), sigma1*ones(n,1), nt, CorrMat);
    prior_runoff_ens(:,:,mm) = w'.*truth.total_runoff;
end

%% Corrupt discharge with error to make measurements

mean1 = 1;
Cv = (0.15)^2; % Cv is a variance
mu1 = log((mean1^2)/sqrt(Cv+mean1^2));
sigma1 = sqrt(log(Cv/(mean1^2)+1));
v = lognrnd(mu1, sigma1, 1e5, 1);
% figure
% histogram(v)
mean(v)
std(v)

error_corrupted_discharge_meas = true_discharge.*lognrnd(mu1, sigma1, nt, m);

g=1; % gage number
ms=10;
figure
plot(tv, error_corrupted_discharge_meas(:,g), 'cyan.', 'markersize', ms)
hold on
plot(tv, true_discharge(:,g), 'k')
legend('measurement','true')
xlabel('Time')
ylabel('Discharge (mm/hr)')

%% Do ensemble ISR

rng(704753262, 'twister')

basin.mask(isnan(basin.mask)) = 0;
basin.mask = logical(basin.mask);
tic
Cv = (0.15)^2; % this is a variance
options.loc = 0;
options.locL = 0.001; % cells
options.locT = 0.001; % hr
options.rho_yy=0;
basin.distmat = calc_distance_matrix(basin.grid_lon(basin.mask), basin.grid_lat(basin.mask));
basin.true_runoff = truth.total_runoff;
basin.t = tv;
basin.true_discharge = true_discharge;
[post_runoff, w1, w2] = ISR_EnKF(prior_runoff_ens, HH, error_corrupted_discharge_meas, s, basin, Cv, options);
toc

mean_prior_runoff = mean(prior_runoff_ens,3);       
sd_prior_runoff = std(prior_runoff_ens, [], 3);
mean_posterior_runoff = mean(post_runoff, 3);
sd_posterior_runoff = std(post_runoff, [], 3);

ens.mean_prior_runoff = mean_prior_runoff;
ens.sd_prior_runoff = sd_prior_runoff;
ens.mean_posterior_runoff = mean_posterior_runoff;
ens.sd_posterior_runoff = sd_posterior_runoff;

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
% ylim([-10,10])

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

for g=1:5
    
    figure
    h1 = plot(tv, mean_prior_discharge(:,g), 'b', 'linewidth', lw);
    hold on
    plot(tv, mean_prior_discharge(:,g) - 2*prior_discharge_sd(:,g), 'b--', 'linewidth', lw);
    plot(tv, mean_prior_discharge(:,g) + 2*prior_discharge_sd(:,g), 'b--', 'linewidth', lw);
    
    h2 = plot(tv, mean_posterior_discharge(:,g), 'r-+', 'linewidth', lw, 'markersize', 10);
    plot(tv, mean_posterior_discharge(:,g) - 2*posterior_discharge_sd(:,g), 'r--', 'linewidth', lw);
    plot(tv, mean_posterior_discharge(:,g) + 2*posterior_discharge_sd(:,g), 'r--', 'linewidth', lw);
    
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

%% Maps (ensemble means)

gi=1:nt;

% basin.mask = flipud(basin.mask);
basin.true_runoff = truth.total_runoff;
plt_ISR_results_overview(basin, mean_prior_runoff, mean_posterior_runoff, truth, tv, gi)

figure
subplot(1,2,1)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, mean_prior_runoff', truth.total_runoff, gi);
title('Prior mean (NSE)')
subplot(1,2,2)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, mean_posterior_runoff', truth.total_runoff, gi);
title('Posterior mean (NSE)')

min(post_runoff(:))

figure
for c = 1:n
    subplot(4,5,c)
    plot(tv, mean_posterior_runoff(c,:),'r')
    hold on
    plot(tv, truth.total_runoff(c,:), 'k')
    if c==1
        legend('posterior','truth')
    end
    title(['Runoff for cell ' num2str(c)])
%     xlim([tv(i1),tv(i2)])
end

%% plot runoff for cells [2,7,12,17] to show benefit of ensemble vs. PW13

% do for time period from May 13-17
i1=378;
i2=i1+24;

% ensemble
figure
ii=1;
for c=[2,7,12,17]
    
    subplot(2,4,ii)
    h1 = plot(tv(i1:i2), ens.mean_posterior_runoff(c,i1:i2),'r', 'linewidth', lw);
    hold on
    plot(tv(i1:i2), ens.mean_posterior_runoff(c,i1:i2) + ens.sd_posterior_runoff(c,i1:i2),'r--', 'linewidth', lw);
    plot(tv(i1:i2), ens.mean_posterior_runoff(c,i1:i2) - ens.sd_posterior_runoff(c,i1:i2),'r--', 'linewidth', lw);
    
    h2 = plot(tv(i1:i2), ens.mean_prior_runoff(c,i1:i2), 'b', 'linewidth', lw);
    plot(tv(i1:i2), ens.mean_prior_runoff(c,i1:i2) + ens.sd_prior_runoff(c,i1:i2), 'b--', 'linewidth', lw);
    plot(tv(i1:i2), ens.mean_prior_runoff(c,i1:i2) - ens.sd_prior_runoff(c,i1:i2), 'b--', 'linewidth', lw);
    h3 = plot(tv(i1:i2), truth.total_runoff(c,i1:i2), 'k', 'linewidth', lw);
    if c==1
        legend([h2,h1,h3],{'Prior','Posterior','Truth'})
    end
    xlim([tv(i1),tv(i2)])
    title(['Runoff for cell ' num2str(c)])
    
    ii=ii+1;
    
end

% PW13
ii=1;
for c=[2,7,12,17]
    
    subplot(2,4,ii+4)
    h1 = plot(tv(i1:i2), pw13.post_runoff(c,i1:i2),'r', 'linewidth', lw);
    hold on
    h2 = plot(tv(i1:i2), pw13.prior_runoff(c,i1:i2), 'b', 'linewidth', lw);
    h3 = plot(tv(i1:i2), truth.total_runoff(c,i1:i2), 'k', 'linewidth', lw);
    if c==1
        legend([h2,h1,h3],{'Prior','Posterior','Truth'})
    end
    xlim([tv(i1),tv(i2)])
    title(['Runoff for cell ' num2str(c)])
    
    ii=ii+1;
    
end

%% Maps for specific cells, different values of alpha

figure(3)
ii=1; % 1, 4, 7
i1 = find(tv==datetime(2006,5,1));
i2 = find(tv==datetime(2006,5,6));
for c=[2,11,17]
    subplot(1,3,ii)
    h1 = plot(tv(i1:i2), mean_prior_runoff(c,i1:i2), 'b','linewidth',lw);
    hold on
    plot(tv(i1:i2), mean_prior_runoff(c,i1:i2) + 2*sd_prior_runoff(c,i1:i2), 'b--','linewidth',lw);
    plot(tv(i1:i2), mean_prior_runoff(c,i1:i2) - 2*sd_prior_runoff(c,i1:i2), 'b--','linewidth',lw);
    h2 = plot(tv(i1:i2), mean_posterior_runoff(c,i1:i2),'r','linewidth',lw);
    plot(tv(i1:i2), mean_posterior_runoff(c,i1:i2) + 2*sd_posterior_runoff(c,i1:i2),'r--','linewidth',lw);
    plot(tv(i1:i2), mean_posterior_runoff(c,i1:i2) - 2*sd_posterior_runoff(c,i1:i2),'r--','linewidth',lw);
    h3 = plot(tv(i1:i2), truth.total_runoff(c,i1:i2), 'k','linewidth',lw);
    ii = ii+1;
    title(['Cell ' num2str(c) ', \alpha = ' num2str(alpha1)])
    xlabel('Time')
    ylabel('Runoff (mm/hr)')
    switch c
        case 2
            ylim([-1, 4]);
            if ii==1, legend('prior','posterior','truth'), end
        case 11
            ylim([-1, 6]);
        case 17
            ylim([-0.2, 1.7]);
    end
    legend([h1,h2,h3],{'Prior','Posterior','Truth'})
    set(gca, 'fontsize', fs)
end

%%
myKGE(truth.total_runoff(c,i1:i2)', mean_posterior_runoff(c,i1:i2)')
myKGE(truth.total_runoff(c,i1:i2)', mean_prior_runoff(c,i1:i2)')

i1 = find(tv==datetime(2006,5,1));
i2 = find(tv==datetime(2006,5,6));
figure(4)
g=1; % gage at the outlet
subplot(2,1,2)
h1 = plot(tv(i1:i2), mean_prior_discharge(i1:i2,g), 'b', 'linewidth', lw);
hold on   
plot(tv(i1:i2), mean_prior_discharge(i1:i2,g) + 2*prior_discharge_sd(i1:i2,g), 'b--', 'linewidth', lw);
plot(tv(i1:i2), mean_prior_discharge(i1:i2,g) - 2*prior_discharge_sd(i1:i2,g), 'b--', 'linewidth', lw);
h2 = plot(tv(i1:i2), mean_posterior_discharge(i1:i2,g), 'r-', 'linewidth', lw, 'markersize', ms);
plot(tv(i1:i2), mean_posterior_discharge(i1:i2,g) + 2*posterior_discharge_sd(i1:i2,g), 'r--', 'linewidth', lw);
plot(tv(i1:i2), mean_posterior_discharge(i1:i2,g) - 2*posterior_discharge_sd(i1:i2,g), 'r--', 'linewidth', lw);
h3 = plot(tv(i1:i2), true_discharge(i1:i2,g), 'k', 'linewidth', lw);
h4 = plot(tv(i1:i2), error_corrupted_discharge_meas(i1:i2,g), 'green.', 'markersize', ms);
legend([h1,h2,h3,h4], {'Prior','Posterior','Truth','Meas.'})
title(['\alpha = ' num2str(alpha1)])
xlabel('Time')
ylabel('Discharge(mm/hr)')
set(gca, 'fontsize', fs)

