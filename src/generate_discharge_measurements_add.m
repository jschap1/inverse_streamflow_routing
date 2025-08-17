% Generates discharge measurements, corrupted with error as necessary
%
% 7/27/2023 JRS

function gage_w_error = generate_discharge_measurements_add(true_discharge, mu1, sigma1)

% gage = true_discharge; % should corrupt with error

gage_all = true_discharge;
gage = true_discharge; % should corrupt with error
% gage = true_discharge_w_swot_sampling; % should corrupt with error

% Y20: additive Gaussian error with mu, sigma
[nt, m] = size(gage);

sigma1 = sigma1*gage; % stddev of additive error (15% of truth)
add_err = mu1 + sigma1.*randn(nt,m);
gage_w_error = gage + add_err;

err = true_discharge(:) - gage_w_error(:);
nanmean(err)
nanstd(err)

% figure
% histogram(err)
% why does this histogram not look gaussian?
% because each outcome is from a different distribution, with a different
% standard deviation (sigma varies for each observation)
% If we repeated this many times, then there would be Gaussian dists

% y = rand(1000,50);
% yobs = y + randn(1000, 50);
% figure
% histogram(yobs-y);

return