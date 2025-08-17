% Function version of "tmpa_vs_nldas_errors"
%
% Checks whether simulated errors are as expected
% 9/25/2023 JRS
%
% Used for checking simulated errors for EnsISR
% Needs work. Function is not showing expected results.
%
% Inputs:
% True runoff (nt, n)
% Runoff prior (nt, n, M)
% addmult (choose 0 for additive errors, 1 for multiplicative errors)

function analyze_errors_func_M(runoff_prior, true_runoff, addmult)

lw = 2;
fs = 18;

%% Calculate error
if addmult==0
    disp('Additive error')
    err = runoff_prior - true_runoff;
elseif addmult==1
    disp('Multiplicativeerror')
    err = runoff_prior./true_runoff;
end

%% Check Gaussianity

[n,nt,M] = size(err);

% All marginals are Gaussian (this does not imply the errors are jointly
% Gaussian)
err_subset = reshape(err(10,10,:), M,1);
% err_subset = reshape(err(:,1,:), n*M,1);
% err_subset = reshape(err(1,:,:), nt*M,1);

figure
histogram(err_subset(:), 'normalization', 'pdf')
hold on;
ax = gca;
x = linspace(ax.XLim(1), ax.XLim(2), 1000);
mu = mean(err_subset(:));
sigma = std(err_subset(:));
title('Runoff error for all times at a single cell')
plot(x, normpdf(x, mu, sigma), 'LineWidth', 2)
legend('Runoff error histogram','Normal PDF', 'location', 'best')
ylabel('Prob')
xlabel('Runoff error (mm/day)')
set(gca, 'fontsize', fs)

% tmp = reshape(err(:,1,:), n*M,1);

% figure, qqplot(err_subset)

%% For additive errors, does standard deviation increase with true/prior runoff?

xt = true_runoff; % could use runoff prior or true runoff
% xt = runoff_prior; % could use runoff prior or true runoff

E = [min(xt(:)):0.5:max(xt(:))]; % bin the prior runoff and see what the error dist looks like in each bin
Y = discretize(xt(:),E);
nbin = length(E);
binsd = zeros(nbin-1,1);
binromean = zeros(nbin-1,1);
numobs = zeros(nbin-1,1);
err_reshape = reshape(err, n*nt, M);
runoff_prior_reshape = reshape(runoff_prior, n*nt,M);

for bin=1:nbin-1
    
    binerr = err_reshape(Y==bin,:); % get runoff errors corresponding to each size TMPA runoff
    binmean(bin) = mean(mean(binerr,2));
    binsd(bin) = mean(std(binerr, [], 2));
    numobs(bin) = size(binerr,1);
    
    plotfig=0;
    if plotfig
        figure
        histogram(binerr)
        title(['bin ' num2str(bin)])
        title(['Prior runoff between ' num2str(E(bin)) ' to ' num2str(E(bin+1))])
        xlabel('Runoff error')
        ylabel('Freq'), set(gca, 'fontsize', fs)
    end
    
    % is there a relationship between the variance of the error and the size of
    % the prior runoff?
    binro = runoff_prior_reshape(Y==bin,:);
    binromean(bin) = mean(mean(binro, 2));
end

binromean(numobs<2)=[];
binsd(numobs<2)=[];

figure, plot(binromean, binsd, '.', 'markersize', 20)
xlabel('mean nldas runoff (truth)')
ylabel('error standard deviation')
title('error stddev vs. runoff')
set(gca, 'fontsize', fs)
% proportionality constant relating prior runoff to variance of runoff
% error: fit regression line
c = corr(binromean, binsd);
disp(['Correlation: ' num2str(c)])

% Add linear fit to plot
f = fitlm(table(binromean, binsd));
alpha1 = f.Coefficients.Estimate;

xvec = linspace(min(binromean),max(binromean));
yvec = alpha1(1) + alpha1(2)*xvec;
hold on
plot(xvec, yvec, 'linewidth', lw)
legend('Binned values','Line of best fit')

%% Check autocorrelation structure

% nt=30;
% TC= calc_TC(1/T, 0, 1, nt-1, 36);
% samp = mvnrnd(zeros(nt,1), TC, M)';
% 
% lags = round(nt/10);
% myacf = zeros(lags,nt);
% for mm=1:M
%     myacf(:,mm) = acf(samp(:,mm), lags, 0);
% end
% figure,bar(0:lags, [1;mean(myacf,2)])

lags = round(nt/10);
myacf = zeros(lags,n,M);
for i=1:n
    for mm=1:M
%     figure
        myacf(:,i,mm) = acf(err(i,:,mm)', lags, 0);
    end
    disp(i)
end

% Plot empirical ACF for a grid cell and replicate
figure
bar(myacf(:,1,1))
title('Mean sample ACF for a cell and replicate')

% Plot empirical ACF for a grid cell, averaged over replicates
figure
bar(mean(myacf(:,1,:),3))
title('Mean sample ACF for a cell, averaged over replicates')

myacf_repavg = mean(myacf,3);
myacf_cellavg = mean(myacf_repavg,2);

% Plot empirical and theoretical ACF (for a cell)
T = 5; % days
figure
bar(0:lags, [1; myacf_repavg(:,101)])
hold on
plot(0:lags, exp(-(0:lags)./T), 'r-*', 'linewidth', lw)
xlabel('Time lag (days)')
ylabel('Correlation')
legend('Empirical ACF', 'Theoretical ACF')
title('One cell, T = 5 days')
set(gca, 'fontsize', fs)

% Plot empirical and theoretical ACF (average over all cells)
figure
bar(0:lags, [1; myacf_cellavg])
hold on
plot(0:lags, exp(-(0:lags)./T), 'r-*', 'linewidth', lw)
xlabel('Time lag (days)')
ylabel('Correlation')
legend('Empirical ACF', 'Theoretical ACF')
title('Average over cells, T = 5 days')
set(gca, 'fontsize', fs)

% Variogram comparisons are done separately in R

return

% ACF and variogram
% This is better done in R, especially variogram

% %% Check autocorrelation function
% 
% [nt,n] = size(err);
% lags = round(nt/10);
% myacf = zeros(lags,n);
% for i=1:n
% %     figure
%     myacf(:,i) = acf(err(:,i), lags, 0);
% end
% 
% % Mean ACF over all grid cells
% figure
% subplot(1,2,1)
% plotacf(mean(myacf,2), lags, n)
% title('Mean sample ACF over cells')
% 
% % ACF bar plot with errorbars to represent variability
% subplot(1,2,2)
% boxplot(myacf')
% line([0,lags+1],[0,0], 'Color', 'k', 'linewidth', 1)
% title('Sample ACFs over all grid cells')
% 
% %% Check semivariogram
% 
% % We check variograms in R, since it is better equipped to handle
% % geospatial data with gstat
% 
% return
