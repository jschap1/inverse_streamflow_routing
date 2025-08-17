% Make figures for the ISR paper (with log EnsISR (not used, though))
%
% 11/10/2023 JRS
% This is the latest script for generating figures for the paper

clear
cd('/home/jschap/Dropbox/ISR/inverse_streamflow_routing')
addpath(genpath('./src/'))
load('./allegheny_data/setup/setup-swot-gage.mat');
i1=60; i2=179;
nt=i2-i1+1;
lw=2;
fs=16;
ms=30;

truth = struct();
true_runoff = true_runoff(:,i1:i2);
truth.true_runoff = true_runoff;
truth.total_runoff = true_runoff;
basin.true_runoff = true_runoff;

colors = [0,0,255; % dark blue (Ens prior)
           0, 255, 0; % green (Y20 posterior)
           255,0,0;  % red (Ens posterior)
           0,255,255]/255; % cyan (Y20 prior)
       
tv = 1:120;
gi = (k+1):nt-(k+1);

% low prior (0.5*nldas) (0.25*tmpa)
y20_AG = load('./allegheny_data/results/outlet_only/Y20AG_m0a1L40T5_daily_5_low_unifprior.mat');

% medium prior (tmpa)
% y20_AG = load('./allegheny_data/results/outlet_only/Y20AG_m0a1L40T5_daily_5_med_unifprior.mat');

% high prior (2*nldas) (4*tmpa)
% y20_AG = load('./allegheny_data/results/outlet_only/Y20AG_m0a1L40T5_daily_5_high_unifprior.mat');

mm=1;
y20_AG.post_runoff = y20_AG.post_runoff(:,:,mm);
y20_AG.post_stddev = y20_AG.post_stddev(:,:,mm);
y20_AG.prior_stddev = y20_AG.prior_stddev(:,:,mm);
y20_AG.runoff_prior = y20_AG.runoff_prior(:,:,mm);

%% Look at overview of results

[kge, rmse, nse] = plt_ISR_results_overview(basin, y20_AG.runoff_prior', ...
    y20_AG.post_runoff', truth, tv, gi)
title('Y20 A/G')

%% Plot ensemble of runoff (Y20_AG) with shadederrorbar

kk=40;

figure(99)
subplot(1,3,1)

meas_times = find(~isnan(y20_AG.gage_w_error(gi,2))); 

% Show median/mean plus or minus IQR

errbars_prior = 1.96*y20_AG.prior_stddev(gi,kk); % 
errbars_post = 1.96*y20_AG.post_stddev(gi,kk); % 

h1 = shadedErrorBar(tv(gi),y20_AG.runoff_prior(gi,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),y20_AG.post_runoff(gi,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

nse1 = myNSE(true_runoff(kk,gi)',y20_AG.runoff_prior(gi,kk))
nse2 = myNSE(true_runoff(kk,gi)',y20_AG.post_runoff(gi,kk))

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Low prior'])
% title(['Runoff at cell ' num2str(kk)])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Runoff (mm/day)')

% ylim([-20,30])
% ylim([-1,7]) % medium
% ylim([-4,11]) % low
grid on
set(gca, 'fontsize', fs)
% shows 95% confidence bounds

