%% PW13 ISR
%
% 7/20/2023 JRS
% Also implementing domain localization

clear, clc, close all
% cd /hdd/ISR/inverse_streamflow_routing
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))
addpath(genpath('/Users/jschap/Documents/MATLAB/cbrewer2'))

load('./allegheny_data/setup-2-gage.mat');
load('./allegheny_data/alleg_nldas_nanfilled.mat')
load('./allegheny_data/alleg_tmpa_nanfilled.mat')

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

true_discharge = gage;
[nt,m] = size(true_discharge);

%% Plot tmpa prior vs. nldas (true) runoff (Y20 Figure 6)

% TMPA prior, NLDAS truth

figure
t=75:77;
cmax = 2.5;
for i=1:3
    subplot(2,3,i)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['TMPA runoff (day ' num2str(t(i)) ')'])
    caxis([0,cmax])
    subplot(2,3,3+i)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['NLDAS runoff (day ' num2str(t(i)) ')'])
    caxis([0,cmax])
end
colormap cool

%% Assemble runoff prior (PW13)

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

%% Generate discharge "measurements" (PW13)

% gage = true_discharge; % should corrupt with error

gage_all = true_discharge;
gage = true_discharge; % should corrupt with error
% gage = true_discharge_w_swot_sampling; % should corrupt with error

% Y20: additive Gaussian error with mu, sigma

mu1 = 0; % mean of additive error
sigma1 = 0.00*gage; % stddev of additive error (15% of truth)
add_err = mu1 + sigma1.*randn(nt,m);
gage_w_error = gage + add_err;

% Plot error for several gages
figure
lw = 2;
ind = 1;
for gg=[1,2]
    subplot(1,2,ind)
    plot(tv, gage_all(:,gg), 'linewidth', lw)
    hold on
    plot(tv, gage_w_error(:,gg), 'red.', 'markersize', 20)
    xlabel('Time')
    legend('Actual discharge','Measured discharge')
    ylabel('Discharge (mmd/day)')
    ind = ind + 1;
end

err = true_discharge(:) - gage_w_error(:);
nanmean(err)
nanstd(err)

figure
histogram(err)
% why does this histogram not look gaussian?
% because each outcome is from a different distribution, with a different
% standard deviation (sigma varies for each observation)
% If we repeated this many times, then there would be Gaussian dists

% y = rand(1000,50);
% yobs = y + randn(1000, 50);
% figure
% histogram(yobs-y);

%% Plot discharge "measurements" (PW13)

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

%% Do ISR (PW13)

% runoff_init = ones(nt,n);

s = k+1;
% s = 2*(k+1)-1; % for window of length 32 days
cov = 1; % coefficient of variation
% R = (0.15^2); % meas error covariance
R = 0;

totmean = mean(nldas.runoff_mean);
totsd = std(nldas.runoff_mean);
runoff_init = ones(nt, n)*totmean;

tic
% [post_runoff_both_gages] = ISR_PW13(tmpa_runoff_prior, HH, gage_w_error, s, 'proportional', cov, R);
[post_runoff_ds_gage] = ISR_PW13(tmpa_runoff_prior, HH(1,:,:), gage_w_error(:,1), s, 'proportional', cov, R);
% [post_runoff_us_gage] = ISR_PW13(tmpa_runoff_prior, HH(2,:,:), gage_w_error(:,2), s, 'proportional', cov, R);
% [post_runoff_PW13] = ISR_PW13(runoff_init, HH, gage, s, 'const_diag', totsd, R);
toc

% post_runoff_both_gages
% post_runoff_ds_gage
% post_runoff_us_gage

% save('./allegheny_data/ISR_results_PW13.mat', 'post_runoff_PW13', 's', 'cov', 'R', 'gage')

%% Investigate negative runoff values

% In the ds gage and both gages cases, there are a few (~10) negative
% runoff values. Where do these occur and why?

[aa,bb] = find(post_runoff_ds_gage<0)

% For ds gage case, these happen at cells [75, 94, 112, 113, 150, 151, 166,
% 169] and at times [169, 169, 169, 169, 222, 222, 43, 222].

% Something is special about times 169 and 222 at cells [75, 94, 112, 113]
% and [150, 151, 169], respectively.

% Investigate the update step within ISR

figure
plot(post_runoff_ds_gage)

negsum = 0;
for j=1:n
    negind = find(post_runoff_ds_gage(:,j)<0);
    nneg = length(negind);
    negsum = sum(post_runoff_ds_gage(negind,j));
    if negsum<biggestneg
        n_most_neg = j;
        biggestneg = negsum;
    end
end

figure
plot(post_runoff_ds_gage(:,n_most_neg))


%% Plot the map for the time with negative runoff (169)

t_ind = 169;

runoff_map_w_neg = make_map(basin, post_runoff_ds_gage(t_ind,:));
tr_map = make_map(basin, true_runoff(:,t_ind)');
prior_map = make_map(basin, tmpa_runoff_prior(t_ind,:));
tv(t_ind)

min(post_runoff_ds_gage(t_ind,:))

figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, prior_map, 'Prior runoff (mm/day)')
caxis([-5,20])
subplot(1,3,2)
plotraster(basin.lonv, basin.latv, runoff_map_w_neg, 'Posterior runoff (mm/day)')
caxis([-5,20])
subplot(1,3,3)
plotraster(basin.lonv, basin.latv, tr_map, 'True runoff (mm/day)')
caxis([-5,20])
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
colormap(colorMap)

% figure
% plot(true_runoff(:,222)', post_runoff_both_gages(222,:), 'b.')
% 
% figure
% plot(true_runoff(:,222)', tmpa_runoff_prior(222,:), 'r.')
% 
% corr(true_runoff(:,222), tmpa_runoff_prior(222,:)')
% corr(true_runoff(:,222), post_runoff_both_gages(222,:)')

%% Plot overview of results

truth = struct();
truth.total_runoff = nldas_runoff_true';
truth.true_runoff = nldas_runoff_true';
basin.true_runoff = truth.true_runoff;

gi = (k+1):nt-(k+1);
plt_ISR_results_overview(basin, tmpa_runoff_prior', post_runoff_ds_gage', truth, tv, gi)

%% Compare subbasin cases

% Plot change in NSE from prior to posterior for each case

% basin.mask = flipud(basin.mask);
[nse1, kge, rmse, nsemap_us] = plot_gofmaps(basin, post_runoff_us_gage, truth.total_runoff, gi);
[nse1, kge, rmse, nsemap_ds] = plot_gofmaps(basin, post_runoff_ds_gage, truth.total_runoff, gi);
[nse1, kge, rmse, nsemap_both] = plot_gofmaps(basin, post_runoff_both_gages, truth.total_runoff, gi);
[nse1, kge, rmse, nsemap_prior] = plot_gofmaps(basin, tmpa_runoff_prior, truth.total_runoff, gi);

% plot outline of subdomain (export raster, convert to vector, plot vector)
R1 = makerefmat(min(basin.lonv), min(basin.latv), 1/8,1/8);
geotiffwrite('./allegheny_data/subdomain1_mask.tif', basin.mask_gage(:,:,1), R1)
geotiffwrite('./allegheny_data/subdomain2_mask.tif', basin.mask_gage(:,:,2), R1)
subdomains_vect = shaperead('./allegheny_data/allegheny_subdomains/subdomains.shp');

%% Figure for paper (no prior correlation case)

fmtstr = 'blue-';
col1 = 'white';
lw2 = 4;
figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, flipud(nsemap_us - nsemap_prior), 'Upstream')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
text(-79.75,41.25,'\DeltaNSE = 0', 'color', col1, 'fontsize', fs)
text(-79.5,42,'\DeltaNSE = 0.44', 'color', col1, 'fontsize', fs)
caxis([-1,1])
title('Upstream gage only')
subplot(1,3,2)
plotraster(basin.lonv, basin.latv, flipud(nsemap_ds - nsemap_prior), 'Upstream')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
text(-79.75,41.25,'\DeltaNSE = 0.24', 'color', col1, 'fontsize', fs)
text(-79.5,42,'\DeltaNSE = 0.26', 'color', col1, 'fontsize', fs)
caxis([-1,1])
title('Downstream gage only')
subplot(1,3,3)
plotraster(basin.lonv, basin.latv, flipud(nsemap_both - nsemap_prior), 'Upstream')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
text(-79.75,41.25,'\DeltaNSE = 0.33', 'color', col1, 'fontsize', fs)
text(-79.5,42,'\DeltaNSE = 0.43', 'color', col1, 'fontsize', fs)
caxis([-1,1])
title('Both gages')
colormap(colorMap)

%%
figure
subplot(1,2,1)
plotraster(basin.lonv, basin.latv, flipud(nsemap_both - nsemap_us), 'Both - US')
caxis([-1,1])
title('Both gages - upstream gage')

subplot(1,2,2)
plotraster(basin.lonv, basin.latv, flipud(nsemap_both - nsemap_ds), 'Both - DS')
caxis([-1,1])
title('Both gages - downstream gage')
colormap(colorMap)

greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
colormap(colorMap)
% colormap(bluewhitered)

% zeros should be clearly indicated
% negatives should get a special color

% nn = 21; % must be odd
% cmap = cbrewer2('seq','YlGnBu',nn); % blues at bottom
% colormap([1,0,0;cmap]); % black for negative values

%% Plot runoff time series for grid cells of interest

% cells of interest
kk1 = 60;
kk2 = 205; 

t1 = 90;
t2 = 180;

kk=kk1;
prior_nse = myNSE(truth.total_runoff(kk,t1:t2)', tmpa_runoff_prior(t1:t2,kk));
both_nse = myNSE(truth.total_runoff(kk,t1:t2)', post_runoff_both_gages(t1:t2,kk));
us_nse = myNSE(truth.total_runoff(kk,t1:t2)', post_runoff_us_gage(t1:t2,kk));
ds_nse = myNSE(truth.total_runoff(kk,t1:t2)', post_runoff_ds_gage(t1:t2,kk));

figure
subplot(1,2,1)
plot(tv(t1:t2), tmpa_runoff_prior(t1:t2,kk1), '-r', 'linewidth', lw)
hold on
plot(tv(t1:t2), post_runoff_us_gage(t1:t2,kk1), '--b', 'linewidth', lw)
plot(tv(t1:t2), post_runoff_ds_gage(t1:t2,kk1), '--g', 'linewidth', lw)
plot(tv(t1:t2), post_runoff_both_gages(t1:t2,kk1), '--cyan', 'linewidth', lw)
plot(tv(t1:t2), truth.total_runoff(kk1,t1:t2), '-k', 'linewidth', lw)
legend(['Prior, nse = ' num2str(prior_nse)],['Posterior (us), nse = ' num2str(us_nse)],...
    ['Posterior (ds), nse = ' num2str(ds_nse)],['Posterior (both), nse = ' num2str(both_nse)],'Truth')
title('Cell 1')
ylabel('Runoff (mm/day)')
xlabel('Date')
set(gca, 'fontsize', fs)

kk=kk2;
prior_nse = myNSE(truth.total_runoff(kk,t1:t2)', tmpa_runoff_prior(t1:t2,kk));
both_nse = myNSE(truth.total_runoff(kk,t1:t2)', post_runoff_both_gages(t1:t2,kk));
us_nse = myNSE(truth.total_runoff(kk,t1:t2)', post_runoff_us_gage(t1:t2,kk));
ds_nse = myNSE(truth.total_runoff(kk,t1:t2)', post_runoff_ds_gage(t1:t2,kk));

subplot(1,2,2)
plot(tv(t1:t2), tmpa_runoff_prior(t1:t2,kk2), '-r', 'linewidth', lw)
hold on
plot(tv(t1:t2), post_runoff_us_gage(t1:t2,kk2), '--b', 'linewidth', lw)
plot(tv(t1:t2), post_runoff_ds_gage(t1:t2,kk2), '--g', 'linewidth', lw)
plot(tv(t1:t2), post_runoff_both_gages(t1:t2,kk2), '--cyan', 'linewidth', lw)
plot(tv(t1:t2), truth.total_runoff(kk2,t1:t2), '-k', 'linewidth', lw)
legend(['Prior, nse = ' num2str(prior_nse)],['Posterior (us), nse = ' num2str(us_nse)],...
    ['Posterior (ds), nse = ' num2str(ds_nse)],['Posterior (both), nse = ' num2str(both_nse)],'Truth')
title('Cell 2')
ylabel('Runoff (mm/day)')
xlabel('Date')
set(gca, 'fontsize', fs)

%% Repeat, but for subbasin averages

basinmask = (basin.mask);
basinmask(isnan(basinmask)) = 0;
basinmask = logical(basinmask);
lon = basin.grid_lon(basinmask);
lat = basin.grid_lat(basinmask);

%% Subdomains map

figure
subplot(1,2,1)
subdomain1 = basin.mask_gage(:,:,1) - basin.mask_gage(:,:,2);
plotraster(basin.lonv, basin.latv, subdomain1, 'Sub Domain 1')
hold on 
plot(basin.gage_lon(1), basin.gage_lat(1), '.b', 'markersize', ms)
text(basin.gage_lon(1), basin.gage_lat(1)+1/16, 'Gage 1', 'fontsize', fs)
plot(basin.gage_lon(2), basin.gage_lat(2), '.b', 'markersize', ms)
text(basin.gage_lon(2), basin.gage_lat(2)+1/16, 'Gage 2', 'fontsize', fs)
% plot(lon(kk1)-1/16, lat(kk1)-1/16, '.r', 'markersize', ms)
% text(lon(kk1), lat(kk1), 'Cell 1', 'fontsize', fs)
% plot(lon(kk2)-1/16, lat(kk2)-1/16, '.r', 'markersize', ms)
% text(lon(kk2), lat(kk2), 'Cell 2', 'fontsize', fs)
xlabel('Lon')
ylabel('Lat')
colorbar off

subplot(1,2,2)
plotraster(basin.lonv, basin.latv, basin.mask_gage(:,:,2), 'Sub Domain 2')
hold on
plot(basin.gage_lon(1), basin.gage_lat(1), '.b', 'markersize', ms)
text(basin.gage_lon(1), basin.gage_lat(1)+1/16, 'Gage 1', 'fontsize', fs)
plot(basin.gage_lon(2), basin.gage_lat(2), '.b', 'markersize', ms)
text(basin.gage_lon(2), basin.gage_lat(2)+1/16, 'Gage 2', 'fontsize', fs)
% plot(lon(kk1)-1/16, lat(kk1)-1/16, '.r', 'markersize', ms)
% text(lon(kk1), lat(kk1), 'Cell 1', 'fontsize', fs)
% plot(lon(kk2)-1/16, lat(kk2)-1/16, '.r', 'markersize', ms)
% text(lon(kk2), lat(kk2), 'Cell 2', 'fontsize', fs)
xlabel('Lon')
ylabel('Lat')
colorbar off

%%
% Method for getting indices of cells in each subdomain
[aa,bb]=(find(subdomain1==1));
subdomain1_index = sub2ind(size(basin.mask), aa,bb);

[aa,bb]=(find(basin.mask_gage(:,:,2)==1));
subdomain2_index = sub2ind(size(basin.mask), aa,bb);

% get runoff for grid cells in each subdomain
% Reshape TMPA runoff prior to map

sub1runoff = struct();
sub2runoff = struct();

sub1runoff.prior = get_runoff_for_cells_in_subdomain(tmpa_runoff_prior, subdomain1_index, basin);
sub2runoff.prior = get_runoff_for_cells_in_subdomain(tmpa_runoff_prior, subdomain2_index, basin);

sub1runoff.us = get_runoff_for_cells_in_subdomain(post_runoff_us_gage, subdomain1_index, basin);
sub2runoff.us = get_runoff_for_cells_in_subdomain(post_runoff_us_gage, subdomain2_index, basin);

sub1runoff.ds = get_runoff_for_cells_in_subdomain(post_runoff_ds_gage, subdomain1_index, basin);
sub2runoff.ds = get_runoff_for_cells_in_subdomain(post_runoff_ds_gage, subdomain2_index, basin);

sub1runoff.both= get_runoff_for_cells_in_subdomain(post_runoff_both_gages, subdomain1_index, basin);
sub2runoff.both = get_runoff_for_cells_in_subdomain(post_runoff_both_gages, subdomain2_index, basin);

sub1runoff.true= get_runoff_for_cells_in_subdomain(truth.total_runoff', subdomain1_index, basin);
sub2runoff.true = get_runoff_for_cells_in_subdomain(truth.total_runoff', subdomain2_index, basin);

% Instead of individual cells, let's plot time series for the average of
% each basin

% Subdomain 1
figure
plot(tv, mean(sub1runoff.prior,1), '-r', 'linewidth', lw)
hold on
plot(tv, mean(sub1runoff.us,1), '--b', 'linewidth', lw)
plot(tv, mean(sub1runoff.ds,1), '--g', 'linewidth', lw)
plot(tv, mean(sub1runoff.both,1), '--cyan', 'linewidth', lw)
plot(tv, mean(sub1runoff.true,1), '-k', 'linewidth', lw)
legend('Prior','Posterior (us) ','Posterior (ds)','Posterior (both)','Truth')
% legend(['Prior, nse = ' num2str(prior_nse)],['Posterior (us), nse = ' num2str(us_nse)],...
%     ['Posterior (ds), nse = ' num2str(ds_nse)],['Posterior (both), nse = ' num2str(both_nse)],'Truth')
title('Subdomain 1')
ylabel('Runoff (mm/day)')
xlabel('Date')
set(gca, 'fontsize', fs)

% Subdomain 2
figure
plot(tv, mean(sub2runoff.prior,1), '-r', 'linewidth', lw)
hold on
plot(tv, mean(sub2runoff.us,1), '--b', 'linewidth', lw)
plot(tv, mean(sub2runoff.ds,1), '--g', 'linewidth', lw)
plot(tv, mean(sub2runoff.both,1), '--cyan', 'linewidth', lw)
plot(tv, mean(sub2runoff.true,1), '-k', 'linewidth', lw)
legend('Prior','Posterior (us) ','Posterior (ds)','Posterior (both)','Truth')
title('Subdomain 2')
ylabel('Runoff (mm/day)')
xlabel('Date')
set(gca, 'fontsize', fs)

% Calculate average NSE for cells in each subdomain

n1 = 233-50;
n2 = 50;

sub1runoff.nse_prior = zeros(n1,1);
for j=1:n1
    sub1runoff.nse_prior(j) = myNSE(sub1runoff.true(j,gi), sub1runoff.prior(j,gi));
    sub1runoff.nse_us(j) = myNSE(sub1runoff.true(j,gi), sub1runoff.us(j,gi));
    sub1runoff.nse_ds(j) = myNSE(sub1runoff.true(j,gi), sub1runoff.ds(j,gi));
    sub1runoff.nse_both(j) = myNSE(sub1runoff.true(j,gi), sub1runoff.both(j,gi));
end

median(sub1runoff.nse_us') - median(sub1runoff.nse_prior')
median(sub1runoff.nse_ds') - median(sub1runoff.nse_prior')
median(sub1runoff.nse_both') - median(sub1runoff.nse_prior')

sub2runoff.nse_prior = zeros(n2,1);
for j=1:n2
    sub2runoff.nse_prior(j) = myNSE(sub2runoff.true(j,gi), sub2runoff.prior(j,gi));
    sub2runoff.nse_us(j) = myNSE(sub2runoff.true(j,gi), sub2runoff.us(j,gi));
    sub2runoff.nse_ds(j) = myNSE(sub2runoff.true(j,gi), sub2runoff.ds(j,gi));
    sub2runoff.nse_both(j) = myNSE(sub2runoff.true(j,gi), sub2runoff.both(j,gi));
end

median(sub2runoff.nse_us') - median(sub2runoff.nse_prior')
median(sub2runoff.nse_ds') - median(sub2runoff.nse_prior')
median(sub2runoff.nse_both') - median(sub2runoff.nse_prior')