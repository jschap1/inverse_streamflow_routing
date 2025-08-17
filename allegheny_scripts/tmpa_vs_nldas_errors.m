%% What do the actual errors look like?
%
% 12/5/2023 JRS
% Also plots true (NLDAS) and prior mean (TMPA) runoff and discharge for
% ISR paper

clear, clc, close all
% cd /Volumes/HD3/ISR/inverse_streamflow_routing
cd /home/jschap/Dropbox/ISR/inverse_streamflow_routing
addpath(genpath('/home/jschap/Documents/MATLAB'))
addpath(genpath('./src/'))

fs=18;lw=2;ms=15;
load('./allegheny_data/setup/setup-swot-gage.mat');
aa=load('./allegheny_data/setup/setup-2-gage.mat');

B = load('./allegheny_data/setup/setup-233-gage.mat');
basin = B.basin;
gage = B.gage;
HH = B.HH;

basin.distmat = aa.basin.distmat;
load('./allegheny_data/setup/alleg_nldas_nanfilled.mat')
load('./allegheny_data/setup/alleg_tmpa_nanfilled.mat')

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

clearvars aa B

%% Cut down to a shorter time period

i1=60; i2=179;
tv = tv(i1:i2);
gage = gage(i1:i2,:);
nldas_runoff = nldas_runoff(:,:,i1:i2);
true_runoff = true_runoff(:,i1:i2); 
[nt,m] = size(gage);
nldas.runoff  = nldas.runoff(:, :, i1:i2);
tmpa.runoff  = tmpa.runoff(:, :, i1:i2);

%% Assemble runoff prior (Y20)

figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be right side up

basin_mask_linear = basin.mask(:);
basin_mask_linear(isnan(basin_mask_linear)) = 0;

tmpa_runoff_linear = reshape(tmpa.runoff, length(basin_mask_linear), nt);
tmpa_runoff_prior = tmpa_runoff_linear(logical(basin_mask_linear),:)';
tmpa_runoff_prior(isnan(tmpa_runoff_prior)) = 0;

nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_true = nldas_runoff_linear(logical(basin_mask_linear),:)';
nldas_runoff_true(isnan(nldas_runoff_true)) = 0;

%% Calculate true and prior discharge

% Calculate it at "river pixels" above flowacc threshold (must include us
% gauge)
% Calculate it at the outlet and the upstream gauge (29, 146)

res = 1/8;
gi = (k+1):nt-(k+1);

figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum)/1e6, 'facc')
hold on
plot(basin.gage_lon(146)-res/2, basin.gage_lat(146)-res/2, '.r','markersize', 20)

% Get grid cells on the river channel above a threshold flowacc
facc_linear = reshape(flipud(basin.flow_accum), length(basin_mask_linear), 1);
facc = facc_linear(logical(basin_mask_linear),:)';
facc(isnan(facc)) = 0;

thres = 5000e6;
ind = facc(:)>thres;
size(HH(ind,:,:))

% Check
facc_map = make_map(basin, facc);
figure
plotraster(basin.lonv, basin.latv, facc_map/1e6, 'facc')

% Calculate true and prior discharge at each river pixel
rivgage.prior_discharge = state_model_dumb(tmpa_runoff_prior, HH(ind, :, :));
rivgage.true_discharge = state_model_dumb(nldas_runoff_true, HH(ind, :, :));
% indices of outlet (29/233) and upstream gauge (146/233): ind(29),
% ind(146) -- can use these to get the specific Q of interest

% Calculate true and prior discharge at the outlet and upstream gage
twogage.prior_discharge = state_model_dumb(tmpa_runoff_prior, HH([29;146], :, :));
twogage.true_discharge = state_model_dumb(nldas_runoff_true, HH([29;146], :, :));

% Convert mm/day to cms for figures
basin_area = max(basin.gage_area)/1000^2; % basin area = 33751 km2
f = (basin_area/n)*1000/86400; % cubic meters per second conversion factor

% Plot true and prior discharge at two gages

twogage.rmse(1) = myRMSE(f*twogage.true_discharge(gi,1), f*twogage.prior_discharge(gi,1));
twogage.rmse(2) = myRMSE(f*twogage.true_discharge(gi,2), f*twogage.prior_discharge(gi,2));
twogage.nse(1) = myNSE(twogage.true_discharge(gi,1), twogage.prior_discharge(gi,1));
twogage.nse(2) = myNSE(twogage.true_discharge(gi,2), twogage.prior_discharge(gi,2));

% define bias as mean(prior) - mean(truth). 
% if bias is positive, the prior is larger than the truth
twogage.bias = nanmean(twogage.prior_discharge - twogage.true_discharge, 1)'; 

%% True vs. prior discharge overview figure

charlbl =  compose("(%s)",('a':'z').'); % labels

figure
subplot(2,2,4)
plot(tv, f*twogage.prior_discharge(:,2), 'b-', 'linewidth', lw)
hold on
plot(tv, f*twogage.true_discharge(:,2), 'k-', 'linewidth', lw)
legend('TMPA prior','NLDAS truth')
title('Upstream gage')
xlabel('Day')
ylabel('Discharge (m^3/s)')
text(0.025,0.95,charlbl{4},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
text(60, 200, ['RMSE = ' num2str(twogage.rmse(2), 3) ' m^3/s'])
text(60, 215, ['NSE = ' num2str(twogage.nse(2), 2)])
text(60, 230, ['Bias = ' num2str(twogage.bias(2), 2) ' m^3/s'])
set(gca, 'fontsize', fs)

subplot(2,2,3)
plot(tv, f*twogage.prior_discharge(:,1), 'b-', 'linewidth', lw)
hold on
plot(tv, f*twogage.true_discharge(:,1), 'k-', 'linewidth', lw)
legend('TMPA prior','NLDAS truth')
title('Outlet')
xlabel('Day')
ylabel('Discharge (m^3/s)')
text(60, 1100, ['RMSE = ' num2str(twogage.rmse(1), 3) ' m^3/s'])
text(60, 1200, ['NSE = ' num2str(twogage.nse(1), 2)])
text(60, 1300, ['Bias = ' num2str(twogage.bias(1), 2) ' m^3/s'])
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
set(gca, 'fontsize', fs)

subplot(2,2,2) % plot of all river pixels' Q with highlighted cells of interest
plot(f*rivgage.true_discharge(:), f*rivgage.prior_discharge(:), '.b', 'markersize', ms)
axis([0,1000,0,1000])
x = linspace(0, 1000, 100);
y = x;
hold on
plot(x, y, 'k-', 'LineWidth', 1);
title('Discharge at 31 river pixels')
grid on
rtd = f*rivgage.true_discharge(gi,:);
rpd = f*rivgage.prior_discharge(gi,:);
rmse = myRMSE(rtd(:), rpd(:));
nse = myNSE(rtd(:), rpd(:));
bias = nanmean(rpd(:) - rtd(:));
%text(600, 100, ['RMSE = ' num2str(rmse, 3) ' m^3/s'])
%text(600, 150, ['NSE = ' num2str(nse, 2)])
%text(600, 200, ['Bias = ' num2str(bias, 2) ' m^3/s'])
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
xlabel('NLDAS truth (m^3/s)')
ylabel('TMPA prior (m^3/s)')
set(gca, 'fontsize', fs)

for i=1:sum(ind)
    nse(i) = myNSE(f*rivgage.true_discharge(:,i), f*rivgage.prior_discharge(:,i));
    rmse(i) = myRMSE(f*rivgage.true_discharge(:,i), f*rivgage.prior_discharge(:,i));
    bias(i) = nanmean(f*rivgage.prior_discharge(:,i) - f*rivgage.true_discharge(:,i));
end
    
% subplot(2,3,4) % histogram of RMSE
% histogram(rmse)
% title('RMSE for 31 river pixels')
% xlabel('RMSE (m^3/s)')
% ylabel('Frequency')
% set(gca, 'fontsize', fs)
% 
% subplot(2,3,5) % histogram of NSE
% histogram(nse)
% title('NSE for 31 river pixels')
% xlabel('NSE')
% ylabel('Frequency')
% set(gca, 'fontsize', fs)

% subplot(2,3,6) % histogram of bias
% histogram(bias)
% title('Bias for 31 river pixels')
% xlabel('Bias (m^3/s)')
% ylabel('Frequency')
% set(gca, 'fontsize', fs)

subplot(2,2,1) % map showing river pixels and highlighting outlet and upstream gage
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum)/1e6, 'Upstream area (km^2)');
hold on
h1 = plot(basin.gage_lon(ind) - res/2, basin.gage_lat(ind) - res/2, 'r.', 'markersize', ms);
h2 = plot(basin.gage_lon(29) - res/2, basin.gage_lat(29) - res/2, 'k.', 'markersize', ms);
h3 = plot(basin.gage_lon(146) - res/2, basin.gage_lat(146) - res/2, 'green.', 'markersize', ms);
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
legend([h2,h3,h1], 'Outlet', 'Upstream gage','River pixels')

%% True vs. prior runoff

% Runoff snapshot
dd = 50; % 10, 50
yupp = max(max(max(nldas.runoff(:,:,dd))), max(max(tmpa.runoff(:,:,dd))));
ylow = min(min(min(nldas.runoff(:,:,dd))), min(min(tmpa.runoff(:,:,dd))));
figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,dd), ['True runoff, Day ' num2str(dd) ' (mm/day)'])
caxis([ylow,yupp])
subplot(1,3,2)
plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,dd), ['Prior runoff, Day ' num2str(dd) ' (mm/day)'])
caxis([ylow,yupp])
ax = subplot(1,3,3);
plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,dd) - nldas.runoff(:,:,dd), 'Prior - Truth (mm/day)')
%colormap cool
% make something similar for runoff, and for runoff error
colormap(ax,bluewhitered)

% Mean average runoff
figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, mean(nldas.runoff,3), 'Mean true runoff (mm/day)')
caxis([0.6,2.3])
subplot(1,3,2)
plotraster(basin.lonv, basin.latv, mean(tmpa.runoff,3), 'Mean prior runoff (mm/day)')
caxis([0.6,2.3])
ax = subplot(1,3,3);
plotraster(basin.lonv, basin.latv, mean(tmpa.runoff,3) - mean(nldas.runoff,3), 'Prior - Truth (mm/day)')
%colormap cool
% make something similar for runoff, and for runoff error
colormap(ax,bluewhitered)

%% If we assume alpha=1, is there anywhere the truth is not within 2sd of the prior?

alpha1 = 1;
tmpa_runoff_sd = alpha1*tmpa.runoff;
ind = struct();
ind.smallprior = nldas.runoff > tmpa.runoff + tmpa_runoff_sd; % prior is too small vs. truth
ind.largeprior = nldas.runoff < tmpa.runoff - tmpa_runoff_sd; % prior is too large vs. truth

% truth is never smaller than 2SD of tmpa_runoff (tmpa_runoff -
% alpha(tmpa_runoff) = 0 for alpha=1
% truth is sometimes bigger than 2SD of tmpa_runoff

sum(ind.smallprior(:))
sum(ind.largeprior(:))

% tmpa runoff is too big sometimes. But it is still guaranteed that truth
% is within 2sd of it

%% Additive error assumption

err = tmpa_runoff_prior - nldas_runoff_true;
figure,histogram(err(:))
xlim([-3,3])
title('TMPA - NLDAS')
xlabel('runoff error (mm/day)')
ylabel('freq')
set(gca, 'fontsize', fs)

% How close are runoff errors to Gaussian?
[f,x,flo,fup] = ecdf(err(:));
ecdf(err(:),'Bounds','on')
mu = mean(err(:));
sigma = std(err(:));
ncdf = normcdf(x,mu,sigma);
figure
plot(x, f, 'linewidth', lw)
hold on
plot(x, ncdf, 'linewidth', lw)
legend('Empirical CDF', 'Normal CDF')
xlabel('x')
ylabel('Pr(x)')
xlim([-3,3])
set(gca, 'fontsize', fs)

E = [0:0.5:5]; % bin the prior runoff and see what the error dist looks like in each bin
Y = discretize(tmpa_runoff_prior(:),E);
nbin = length(E);
binsd = zeros(nbin,1);
binromean = zeros(nbin,1);
for bin=1:nbin-1
    
    binerr = err(Y==bin); % get runoff errors corresponding to each size TMPA runoff
    
%     figure,histogram(binerr), title(['bin ' num2str(bin)])
%     title(['TMPA runoff between ' num2str(E(bin)) ' to ' num2str(E(bin+1))])
%     xlabel('Runoff error')
%     ylabel('Freq'), set(gca, 'fontsize', fs)
    
    binmean(bin) = mean(binerr);
    binsd(bin) = std(binerr);

    % is there a relationship between the variance of the error and the size of
    % the prior runoff?
    binro = tmpa_runoff_prior(Y==bin);
    binromean(bin) = mean(binro);
end

f = fitlm(table(binromean, binsd),'intercept', false);
 % proportionality constant relating 
alpha1 = f.Coefficients.Estimate;

fs=18;
figure, plot(binromean, binsd, '.', 'markersize', 25)
xlabel('TMPA runoff (mm/day)')
ylabel('Error standard deviation (mm/day)')
title('Runoff error stddev vs. TMPA runoff')
text(1.5,0.5,['error stddev = ' num2str(alpha1(end), 2) ' *(TMPA runoff)'], 'fontsize', 14)
grid on
set(gca, 'fontsize', fs)
% proportionality constant relating prior runoff to variance of runoff
% error: fit regression line
% corr(binromean, binsd)

xvec = linspace(min(binromean),max(binromean));
yvec = alpha1(end)*xvec;
% yvec = alpha1(1) + alpha1(end)*xvec;
hold on
ylim([0,2])
plot(xvec, yvec, 'linewidth', lw)

%% Multiplicative error assumption

err = tmpa_runoff_prior./nldas_runoff_true;
figure,histogram(err(:))

E = [0:0.5:5];
Y = discretize(tmpa_runoff_prior(:),E);
nbin = length(E);
binsd = zeros(nbin,1);
binromean = zeros(nbin,1);
for bin=1:nbin-1
    binerr = err(Y==bin);
    figure,histogram(binerr)
    binmean(bin) = mean(binerr);
    binsd(bin) = std(binerr);

    % is there a relationship between the variance of the error and the size of
    % the prior runoff?
    binro = tmpa_runoff_prior(Y==bin);
    binromean(bin) = mean(binro);
end

figure, plot(binromean, binsd, '.', 'markersize', 20)
xlabel('mean tmpa runoff (prior)')
ylabel('error standard deviation')
% proportionality constant relating prior runoff to variance of runoff
% error: fit regression line
corr(binromean, binsd)

f = fitlm(table(binromean, binsd));
 % proportionality constant relating 
alpha1 = f.Coefficients.Estimate(2);