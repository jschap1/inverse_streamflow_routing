% Copied from make_figs. Include as needed in make_figs_E1_4.m

%%

figure
subplot(3,1,1) % Y20 AG
plot(E1.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
E1.kge_post = plot_discharge_ts(E1.postQ(gi,gagei), true_discharge(gi,gagei));
title('Y20 A/G')
plot(E1.gage_w_error(gi,gagei), 'k.', 'markersize', ms)
legend('Prior','Truth','Posterior','Observations')
set(gca, 'fontsize', fs)

subplot(3,1,2) % ENS AG
plot(E4.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
E4.kge_post = plot_discharge_ts(E4.postQ(gi,gagei), true_discharge(gi,gagei));
title('Ens A/G')
plot(E4.error_corrupted_discharge_meas(gi,gagei), 'k.', 'markersize', ms)
set(gca, 'fontsize', fs)
legend('off')

subplot(3,1,3) % ENS LM
plot(E3.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
E3.kge_post = plot_discharge_ts(E3.postQ(gi,gagei), true_discharge(gi,gagei));
title('Ens L/M')
plot(E3.error_corrupted_discharge_meas(gi,gagei), 'k.', 'markersize', ms)
set(gca, 'fontsize', fs)
legend('off')

%% Low vs. high prior analysis

% Are low priors with multiplicative errors problematic?

% gi=35:45;

% calculate NSE for all cells
nse.posterior;

% find cells where TMPA prior is very low/very high
tr = true_runoff(:,gi)';
prior = E1.tmpa_runoff_prior(gi,:);
post = E1.post_runoff(gi,:);

kk=40; % very high prior
kk=101; % very low prior

% focus on timesteps 35:45 becuase there is often a peak in true runoff there

% for kk=1:10
%     figure
%     plot(y20_AG.tmpa_runoff_prior(gi,kk), 'blue', 'linewidth', lw)
%     hold on
%     plot(y20_AG.post_runoff(gi,kk), 'red', 'linewidth', lw)
%     plot(tr(gi,kk), 'k', 'linewidth', lw)
%     legend('Prior','Posterior','Truth')
%     title(['Cell ' num2str(kk)])
% end

% how high or low is the prior compared to the truth?
% use one-sided ER95 to see if truth is contained by the prior
z = 1.645;
figure
plot(prior(:,k) - z*E1.prior_stddev(gi,k))
hold on,plot(tr(:,k),'k')

% how do I know if prior is big? small?
% If prior plus 1.645*sd is smaller than truth, prior is small
% If prior minus 1.645*sd is larger than truth, prior is large
z=1.645;
small_ind = prior + z*E1.prior_stddev(gi,:) < tr;
% large_ind = ~small_ind;
large_ind = prior - z*E1.prior_stddev(gi,:) > tr;
% small_ind = prior < tr;
large_ind = prior > tr;

figure
imagesc(small_ind) % there are a good number of small prior locations
% % but very few places where the prior is large, by this metric
% 
% % how good are the estimates where the prior is small vs. other estimates?
% myNSE(tr(small_ind), y20_AG.runoff_prior(small_ind))
% myNSE(tr(~small_ind), y20_AG.runoff_prior(~small_ind))
% 
% myKGE(tr(small_ind), y20_AG.post_runoff(small_ind))
% myKGE(tr(~small_ind), y20_AG.post_runoff(~small_ind))

%% Revised scatterplots figure

figure
subplot(1,2,1)
plot(tr(~small_ind), prior(~small_ind), '.red')
hold on
plot(tr(small_ind), prior(small_ind), '.blue')
ylabel('Estimate')
title('Prior runoff')
legend('Not small prior','Small prior')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'black', 'LineStyle', '--');

subplot(1,2,2)
plot(tr(~small_ind), post(~small_ind), '.red')
hold on
plot(tr(small_ind), post(small_ind), '.blue')
ylabel('Estimate')
title('Posterior runoff')
legend('Not small prior','Small prior')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'black', 'LineStyle', '--');

notsmallnseprior = myNSE(tr(~small_ind), prior(~small_ind))
notsmallnsepost = myNSE(tr(~small_ind), post(~small_ind))
smallnseprior = myNSE(tr(small_ind), prior(small_ind))
smallnsepost = myNSE(tr(small_ind), post(small_ind))


%% Scatterplots figure (draft)

figure

subplot(2,3,1)
plot(tr(small_ind), prior(small_ind), '.')
ylabel('Estimate')
title('Prior runoff (small prior)')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,2)
plot(tr(large_ind), prior(large_ind), '.')
ylabel('Estimate')
title('Prior runoff (large prior)')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,3)
plot(tr(:), prior(:), '.')
ylabel('Estimate')
title('Prior runoff (all priors)')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,4)
plot(tr(small_ind), post(small_ind), '.')
ylabel('Estimate')
title('Posterior runoff, given small prior')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,5)
plot(tr(large_ind), post(large_ind), '.')
ylabel('Estimate')
title('Posterior runoff, given large prior')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,6)
plot(tr(:), post(:), '.')
ylabel('Estimate')
title('Posterior runoff, given all priors')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

% Color code by size of prior runoff relative to truth
% Can use darker colors for larger difference between prior and truth
% Use same color scheme for prior and posterior scatterplot
% Include NSE values on plots
% Repeat for AGflat method and compare to see effect of proportional errors

% Does this change for the AGflat case? 
% If prior runoff is small, does NSE improve?
% If prior runoff is large, does NSE improve?

% overall
myNSE(tr(:), prior(:))
myNSE(tr(:), post(:))

% not small (similar to overall)
myNSE(tr(~small_ind), prior(~small_ind))
myNSE(tr(~small_ind), post(~small_ind))

myNSE(tr(small_ind), prior(small_ind))
myNSE(tr(small_ind), post(small_ind))

myNSE(tr(large_ind), prior(large_ind))
myNSE(tr(large_ind), post(large_ind))

%% Plot ensemble of runoff (Y20_AG) with shadederrorbar

kk=93; % has negative runoff
kk=44; 
kk=101;
kk=40;

figure
%subplot(1,3,2)

meas_times = find(~isnan(E1.gage_w_error(gi,2))); 

% Show median/mean plus or minus IQR

errbars_prior = 1.96*E1.prior_stddev(gi,kk); % 
errbars_post = 1.96*E1.post_stddev(gi,kk); % 

h1 = shadedErrorBar(tv(gi),E1.runoff_prior(gi,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),E1.post_runoff(gi,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

nse1 = myNSE(true_runoff(kk,gi)',E1.runoff_prior(gi,kk))
nse2 = myNSE(true_runoff(kk,gi)',E1.post_runoff(gi,kk))

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

% plot([meas_times meas_times], [-5,20], 'k--', 'linewidth', 1)

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Medium prior'])
% title(['Runoff at cell ' num2str(kk)])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Runoff (mm/day)')

% ylim([-20,30])
% ylim([-1,7])
% ylim([-4,11])
grid on
set(gca, 'fontsize', fs)
% shows 95% confidence bounds

%% Plot ensemble of runoff (Y20_AG) with lines

z = 1.96;

figure
h1 = plot(E1.runoff_prior(:,kk), 'b*-', 'linewidth', lw);
hold on
plot(E1.runoff_prior(:,kk) + z*E1.prior_stddev(:,kk), 'b-', 'linewidth', lw)
plot(E1.runoff_prior(:,kk) - z*E1.prior_stddev(:,kk), 'b-', 'linewidth', lw)
h2 = plot(E1.post_runoff(:,kk), 'r*-', 'linewidth', lw);
plot(E1.post_runoff(:,kk) + z*E1.post_stddev(:,kk), 'r-', 'linewidth', lw)
plot(E1.post_runoff(:,kk) - z*E1.post_stddev(:,kk), 'r-', 'linewidth', lw)
h3 = plot(true_runoff(kk,:), 'k-', 'linewidth', lw);
legend([h1,h2,h3], 'Prior','Posterior','Truth')
xlabel('t')
ylabel('runoff (mm/day)')

title('Y20, A/G')
set(gca, 'fontsize', fs)

myKGE(true_runoff(kk,:)', E1.runoff_prior(:,kk))
myKGE(true_runoff(kk,:)', E1.post_runoff(:,kk))

myNSE(true_runoff(kk,:)', E1.runoff_prior(:,kk))
myNSE(true_runoff(kk,:)', E1.post_runoff(:,kk))

myRMSE(true_runoff(kk,:)', E1.runoff_prior(:,kk))
myRMSE(true_runoff(kk,:)', E1.post_runoff(:,kk))

%% Plot ensemble of runoff (Y20_LM) with shadederrorbar

kk=40;

figure
%subplot(2,1,1)

meas_times = find(~isnan(E2.gage_w_error(gi,2))); 

% Show median/mean plus or minus IQR

errbars_prior = 1.96*E2.prior_stddev(gi,kk); % 
errbars_post = 1.96*E2.post_stddev(gi,kk); % 

h1 = shadedErrorBar(tv(gi),E2.runoff_prior(gi,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),E2.post_runoff(gi,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

nse1 = myNSE(true_runoff(kk,gi)',E2.runoff_prior(gi,kk))
nse2 = myNSE(true_runoff(kk,gi)',E2.post_runoff(gi,kk))

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

% plot([meas_times meas_times], [-5,20], 'k--', 'linewidth', 1)

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Runoff at cell ' num2str(kk)])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Runoff (mm/day)')

% ylim([-6,20])
grid on
set(gca, 'fontsize', fs)
% shows 95% confidence bounds

%% Plot ensemble of runoff (Ens_AG) with shadederrorbar

% Show median +- IQR
E4.prior_runoff_med = median(E4.prior_runoff_ens,3)';
E4.post_runoff_med = median(E4.post_runoff,3)';
E4.prior_runoff_IQR = iqr(E4.prior_runoff_ens,3)';
E4.post_runoff_IQR = iqr(E4.post_runoff,3)';
E4.prior_runoff_std = std(E4.prior_runoff_ens, [], 3)';
E4.post_runoff_std = std(E4.post_runoff, [], 3)';
% Better yet, show median, 2.5th percentile, and 97.5th percentile
% to reflect UR95 better and avoid making it look like there can be negative values 

E4.prior975 = prctile(E4.prior_runoff_ens, 97.5, 3)';
E4.prior25 = prctile(E4.prior_runoff_ens, 2.5, 3)';
E4.post975 = prctile(E4.post_runoff, 97.5, 3)';
E4.post25 = prctile(E4.post_runoff, 2.5, 3)';

% errbars_prior = [ens_AG.prior975(:,kk), ens_AG.prior25(:,kk)];
% errbars_post = [ens_AG.post975(:,kk), ens_AG.post25(:,kk)];

% errbars_prior = ens_AG.prior_runoff_IQR(:,kk);
% errbars_post = ens_AG.post_runoff_IQR(:,kk);

figure
subplot(2,1,2)
kk=93;
%kk=44;

errbars_prior = 1.96*E4.prior_runoff_std(gi,kk);
errbars_post = 1.96*E4.post_runoff_std(gi,kk);

h1 = shadedErrorBar(tv(gi),E4.meanpriorrunoff(gi,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),E4.mean_post_runoff(gi,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Runoff at cell ' num2str(kk)])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
ylim([-6,20])
grid on
set(gca, 'fontsize', fs)

myKGE(true_runoff(kk,gi)', E4.meanpriorrunoff(gi,kk))
myKGE(true_runoff(kk,gi)', E4.mean_post_runoff(gi,kk))

myNSE(true_runoff(kk,gi)', E4.meanpriorrunoff(gi,kk))
myNSE(true_runoff(kk,gi)', E4.mean_post_runoff(gi,kk))

myRMSE(true_runoff(kk,gi)', E4.meanpriorrunoff(gi,kk))
myRMSE(true_runoff(kk,gi)', E4.mean_post_runoff(gi,kk))

%% GoF statistics

% Change in runoff RMSE for all cells for each method

%% Plot ensemble of runoff (Ens_LM) with shadederrorbar

% Show median +- IQR
E3.prior_runoff_med = median(E3.prior_runoff_ens,3)';
E3.post_runoff_med = median(E3.post_runoff,3)';
E3.prior_runoff_IQR = iqr(E3.prior_runoff_ens,3)';
E3.post_runoff_IQR = iqr(E3.post_runoff,3)';

% Better yet, show median, 2.5th percentile, and 97.5th percentile
% to reflect UR95 better and avoid making it look like there can be negative values 

errbars_prior = E3.prior_runoff_IQR(:,kk);
errbars_post = E3.post_runoff_IQR(:,kk);

figure
h1 = shadedErrorBar(tv,E3.prior_runoff_med(:,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,E3.post_runoff_med(:,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3 = plot(true_runoff(kk,:), 'k-', 'linewidth', lw);

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Median','Posterior Median',...
    'Truth');

title('Ens runoff estimates at grid cell 93')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
set(gca, 'fontsize', fs)

% Quantify uncertainty shrinkage using IQR -> UR95

%% Plot ensemble of runoff (Ens_LM) with lines

figure
h1 = plot(E3.meanpriorrunoff(:,kk), 'b*-', 'linewidth', lw);
hold on
plot(E3.prior_runoff_2p5(:,kk), 'b-', 'linewidth', lw)
plot(E3.prior_runoff_97p5(:,kk), 'b-', 'linewidth', lw)
h2 = plot(E3.mean_post_runoff(:,kk), 'r*-', 'linewidth', lw);
plot(E3.post_runoff_2p5(:,kk), 'r-', 'linewidth', lw)
plot(E3.post_runoff_97p5(:,kk), 'r-', 'linewidth', lw)
h3 = plot(true_runoff(kk,:), 'k-', 'linewidth', lw);
legend([h1,h2,h3], 'Prior','Posterior','Truth')
xlabel('t')
ylabel('runoff (mm/day)')
title('Ens, L/M')
set(gca, 'fontsize', fs)

myKGE(true_runoff(kk,:)', E3.meanpriorrunoff(:,kk))
myKGE(true_runoff(kk,:)', E3.mean_post_runoff(:,kk))

myNSE(true_runoff(kk,:)', E3.meanpriorrunoff(:,kk))
myNSE(true_runoff(kk,:)', E3.mean_post_runoff(:,kk))

myRMSE(true_runoff(kk,:)', E3.meanpriorrunoff(:,kk))
myRMSE(true_runoff(kk,:)', E3.mean_post_runoff(:,kk))

%% Plot discharge with uncertainty (Y20 AG)

gg=2; % outlet

% calculate discharge uncertainty (error propagation)
% y20_AG.priorQ_std = propagate_runoff_error(y20_AG.prior_stddev, HH, gg);
% y20_LM.priorQ_std = propagate_runoff_error(y20_LM.prior_stddev, HH, gg);
% y20_AG.postQ_std = propagate_runoff_error(y20_AG.post_stddev, HH, gg);
% y20_LM.postQ_std = propagate_runoff_error(y20_LM.post_stddev, HH, gg);

E1.priorQ_std = 0;
E2.priorQ_std = 0;
E1.postQ_std = 0;
E2.postQ_std = 0;

figure
h1=plot(E1.priorQ(:,gg), 'b-*', 'linewidth', lw);
hold on
plot(E1.priorQ(:,gg) + E1.priorQ_std, 'b-', 'linewidth', 1);
plot(E1.priorQ(:,gg) - E1.priorQ_std, 'b-', 'linewidth', 1);
h2=plot(E1.postQ(:,gg), 'r-*', 'linewidth', lw);
plot(E1.postQ(:,gg) + E1.postQ_std, 'r-', 'linewidth', 1);
plot(E1.postQ(:,gg) - E1.postQ_std, 'r-', 'linewidth', 1);
h3=plot(true_discharge(:,gg), 'k', 'linewidth', lw);
h4 = plot(E1.gage_w_error(:,gg), 'k.', 'markersize', ms);
legend([h1,h2,h3,h4],'Prior','Posterior','Truth', 'Obs')
xlabel('t')
xlim([20,60])
ylim([0,800])
title('Y20 AG outlet discharge (15% error)')
ylabel('discharge (mm/day)')
set(gca, 'fontsize', fs)

%% Plot discharge with uncertainty (ENS LM)

% It might be worth expressing uncertainty in terms of quantiles instead of
% standard deviation..

% ensAG
[E4.priorQmean, E4.priorQstd, E4.priorQens, E4.priormedianQ,E4.priorq2p5, E4.priorq97p5, E4.priorQiqr] = ...
    calc_ensemble_discharge_moments(E4.prior_runoff_ens, HH);
[E4.postQmean, E4.postQstd, E4.postQens, E4.postmedianQ,E4.postq2p5, E4.postq97p5, E4.postQiqr] = ...
    calc_ensemble_discharge_moments(E4.post_runoff, HH);

% ensLM
[E3.priorQmean, E3.priorQstd, E3.priorQens, E3.priormedianQ,E3.priorq2p5, E3.priorq97p5, E3.priorQiqr] = ...
    calc_ensemble_discharge_moments(E3.prior_runoff_ens, HH);
[E3.postQmean, E3.postQstd, E3.postQens, E3.postmedianQ,E3.postq2p5, E3.postq97p5, E3.postQiqr] = ...
    calc_ensemble_discharge_moments(E3.post_runoff, HH);

%% Comparing distributions of prior and posterior runoff btw Y20 and EnsISR

% Plot PDF for Y20 prior and posterior runoff

runoff_case = 1;
switch runoff_case
    case 1
        tt=110;
        kk=93; % highlighting a time when Y20 posterior runoff<0 and when stddev is reduced significantly
        x = linspace(-50,150);
        ymax = 0.06;
    case 2
        tt=41;
        kk=47; % highlighting a time when Y20 and EnsISR both give similar (accurate) results
        x = linspace(-5,10);
        ymax = 0.6;
end
xmin = min(x);
xmax = max(x);

figure
y = normpdf(x,E1.runoff_prior(tt,kk),E1.prior_stddev(tt,kk));
plot(x,y, 'color', colors(1,:))
hold on
y = normpdf(x,E1.post_runoff(tt,kk),E1.post_stddev(tt,kk));
plot(x,y, 'color', colors(3,:))
plot([true_runoff(kk,tt), true_runoff(kk,tt)], [0, ymax], 'k--', 'linewidth', lw);
title('Y20 runoff estimates')

figure
histogram(E4.prior_runoff_ens(kk,tt,:), 'normalization', 'pdf', 'facecolor', colors(1,:))
hold on
plot([true_runoff(kk,tt), true_runoff(kk,tt)], [0, ymax], 'k--', 'linewidth', lw);
histogram(E4.post_runoff(kk,tt,:), 'normalization', 'pdf', 'facecolor', colors(3,:))
title('EnsISR (with log transform) runoff estimates')
xlim([xmin,xmax])

% Plot histogram for EnsISR prior and posterior runoff

figure
histogram(E3.prior_runoff_ens(kk,tt,:), 'normalization', 'pdf', 'facecolor', colors(1,:))
hold on
plot([true_runoff(kk,tt), true_runoff(kk,tt)], [0, ymax], 'k--', 'linewidth', lw);
histogram(E3.post_runoff(kk,tt,:), 'normalization', 'pdf', 'facecolor', colors(3,:))
title('EnsISR (with log transform) runoff estimates')
xlim([xmin,xmax])

%% Discharge with confidence bounds (Y20)

gg=2;
lw=3;
fs = 22;
z = 1.96;

errbars_prior = zeros(nt,1);
errbars_post = zeros(nt,1);

% errbars_prior = z*y20_AG.prior_stddevQ(:,gg);
% errbars_post = z*y20_AG.post_stddevQ(:,gg);

tvdates = datetime(2009,1,1):datetime(2009,12,31);
tvdates = tvdates(i1:i2);
tv = 1:120;

figure
h1 = shadedErrorBar(tv,E1.priorQ(:,gg),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});

hold on

h2 = shadedErrorBar(tv,E1.postQ(:,gg),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3=plot(E4.priorQmean(:,gg), '-*', 'linewidth', lw, 'color', colors(1,:));
h4=plot(E4.postQmean(:,gg), '-*', 'linewidth', lw, 'color', colors(2,:));

h5 = plot(true_discharge(:,gg), 'k-', 'linewidth', lw);
h6 = plot(E4.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);

% retrace the EnsISR posterior median for visibility
plot(tv,E1.postQ(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));

% xlim([15,85])
legend([h1.mainLine, h2.mainLine,h3,h4,h5,h6], ...
    'Y20 Prior', 'Y20 Posterior', ...
    'Ens Prior','Ens Posterior',...
    'Truth','Observations');
% this plot shows median +- IQR. For normal distribution, median = mean

title('Discharge estimates at Allegheny Basin outlet')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
set(gca, 'fontsize', fs)


%% Discharge with confidence bounds (EnsISR)

gg=2;
lw=3;
fs = 22;
z = 1.96;

errbars_prior = z*E4.priorQstd(:,gg);
errbars_post = z*E4.postQstd(:,gg);

% errbars_prior = [ens_AG.priorq97p5(:,gg),ens_LM.priorq2p5(:,gg)];
% errbars_post = [ens_LM.postq97p5(:,gg),ens_LM.postq2p5(:,gg)];
% errbars_prior = ens_LM.priorQiqr(:,gg);
% errbars_post = ens_LM.postQiqr(:,gg);

tvdates = datetime(2009,1,1):datetime(2009,12,31);
tvdates = tvdates(i1:i2);
tv = 1:120;

figure
h1 = shadedErrorBar(tv,E4.priorQmean(:,gg),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});

hold on

h2 = shadedErrorBar(tv,E4.postQmean(:,gg),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3=plot(E1.priorQ(:,gg), '-*', 'linewidth', lw, 'color', colors(1,:));
h4=plot(E1.postQ(:,gg), '-*', 'linewidth', lw, 'color', colors(2,:));

h5 = plot(true_discharge(:,gg), 'k-', 'linewidth', lw);
h6 = plot(E4.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);

% retrace the EnsISR posterior median for visibility
plot(tv,E4.postQmean(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));

% xlim([15,85])
legend([h1.mainLine, h2.mainLine,h3,h4,h5,h6], ...
    'Ens Prior', 'Ens Posterior', ...
    'Y20 Prior','Y20 Posterior',...
    'Truth','Observations');
% this plot shows median +- IQR. For normal distribution, median = mean

title('Discharge estimates at Allegheny Basin outlet')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
set(gca, 'fontsize', fs)

%% Discharge, with means and standard deviations for EnsISR

figure

h1=plot(E3.priorQmean(:,gg), '-', 'linewidth', lw, 'color', colors(1,:));
hold on
plot(E3.priorQmean(:,gg) + E3.priorQstd(:,gg), '-', 'linewidth', 1, 'color', colors(1,:));
plot(E3.priorQmean(:,gg) - E3.priorQstd(:,gg), '-', 'linewidth', 1, 'color', colors(1,:));

h2=plot(E3.postQmean(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));
plot(E3.postQmean(:,gg) + E3.postQstd(:,gg), '-', 'linewidth', 1, 'color', colors(3,:));
plot(E3.postQmean(:,gg) - E3.postQstd(:,gg), '-', 'linewidth', 1, 'color', colors(3,:));

h5=plot(E1.priorQ(:,gg), '-*', 'linewidth', lw, 'color', colors(1,:));
h6=plot(E1.postQ(:,gg), '-*', 'linewidth', lw, 'color', colors(2,:));

h3=plot(true_discharge(:,gg), 'k-', 'linewidth', lw);
h4 = plot(E3.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);

legend([h5,h6,h1,h2,h3,h4],'Y20 Prior','Y20 Posterior','Ens Prior','Ens Posterior','Truth', 'Obs')
xlabel('t')
% xlim([20,100])
ylim([0,1200])
title('Ens LM outlet discharge (15% error)')
ylabel('discharge (mm/day)')
set(gca, 'fontsize', fs)

%% Show uncertainty ratio for discharge

m=14; % number of gauges
for gg=1:m
    
gi = 4:120; % avoids spin up (no estimate)
top5 = E4.postQmean(gi,gg) + z*E4.postQstd(gi,gg);
bottom5 = E4.postQmean(gi,gg) - z*E4.postQstd(gi,gg);
E4.ur95_post(gg) = sum(top5 - bottom5,1)'./sum(true_discharge(gi,gg));

top5 = E4.priorQmean(gi,gg) + z*E4.priorQstd(gi,gg);
bottom5 = E4.priorQmean(gi,gg) - z*E4.priorQstd(gi,gg);
E4.ur95_prior(gg) = sum(top5 - bottom5,1)'./sum(true_discharge(gi,gg));

end

% Change in UR95 shows how much the 95 percent confidence interval shrunk
delta_UR = E4.ur95_prior - E4.ur95_post
% you want the posterior to be smaller than the prior

figure
boxplot([E4.ur95_prior', E4.ur95_post'])

%% Show 95% exceedence ratio for discharge

E4.er95_prior = myER(true_discharge', permute(E4.postQens,[2,1,3]), 95);
E4.er95_post = myER(true_discharge', permute(E4.priorQens,[2,1,3]), 95);

delta_ER = E4.er95_prior - E4.er95_post
% you want the posterior to be smaller than the prior

figure
boxplot([E4.er95_prior, E4.er95_post])
