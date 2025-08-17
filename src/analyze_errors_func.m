% Function version of "tmpa_vs_nldas_errors"
%
% Checks whether actual errors meet Y20 ISR assumptions
% 9/25/2023 JRS
% Problem is that we only have one realization (TMPA), so we can't really
% check for Gaussianity at a particular place and time. We can check
% assumptions on correlation and the relationship between runoff and runoff
% error variance
%
% Inputs:
% True runoff (nt, n)
% Runoff prior (nt, n, M)
% addmult (choose 0 for additive errors, 1 for multiplicative errors)

function analyze_errors_func(runoff_prior, true_runoff, addmult)

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

% How close are runoff errors to Gaussian? Lognormal?

gamma1 = min(err(:)); % lognormal correction factor
errc = log(err - gamma1 + eps); % avoids negative values, so we can take log

% QQ plots
figure
subplot(2,2,1)
histogram(err(:))
grid on
xlabel('Runoff (mm/day)')
title('Runoff errors')
set(gca, 'fontsize', fs)
subplot(2,2,2)
histogram(errc(:))
xlim([0,3])
grid on
xlabel('log (Runoff (mm/day))')
title('(Shifted) log runoff errors')
set(gca, 'fontsize', fs)
subplot(2,2,3)
qqplot(err(:))
grid on
title('Runoff errors')
set(gca, 'fontsize', fs)
subplot(2,2,4)
qqplot(errc(:))
grid on
ylim([-3,3])
title('Log runoff errors')
set(gca, 'fontsize', fs)

%% For additive errors, does standard deviation increase with true runoff?

xt = true_runoff; % could use runoff prior or true runoff
% xt = runoff_prior; % could use runoff prior or true runoff

E = [min(xt(:)):0.5:max(xt(:))]; % bin the prior runoff and see what the error dist looks like in each bin
Y = discretize(xt(:),E);
nbin = length(E);
binsd = zeros(nbin,1);
binromean = zeros(nbin,1);
numobs = zeros(nbin,1);
for bin=1:nbin-1
    
    binerr = err(Y==bin); % get runoff errors corresponding to each size TMPA runoff
    binmean(bin) = mean(binerr);
    binsd(bin) = std(binerr);
    numobs(bin) = length(binerr);
    
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
    binro = runoff_prior(Y==bin);
    binromean(bin) = mean(binro);
end

figure, plot(binromean, binsd, '.', 'markersize', 20)
xlabel('mean tmpa runoff (prior)')
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

%% Check autocorrelation function

[nt,n] = size(err);
lags = round(nt/10);
myacf = zeros(lags,n);
for i=1:n
%     figure
    myacf(:,i) = acf(err(:,i), lags, 0);
end

% Mean ACF over all grid cells
figure
subplot(1,2,1)
plotacf(mean(myacf,2), lags, n)
title('Mean sample ACF over cells')

% ACF bar plot with errorbars to represent variability
subplot(1,2,2)
boxplot(myacf')
line([0,lags+1],[0,0], 'Color', 'k', 'linewidth', 1)
title('Sample ACFs over all grid cells')

%% Check semivariogram

% We check variograms in R, since it is better equipped to handle
% geospatial data with gstat

return
