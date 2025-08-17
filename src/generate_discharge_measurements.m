% Generate discharge measurements
%
% 6/15/2023 JRS
% Corrupts true discharge with uncorrelated normal or lognormal errors
%
% INPUTS
% m = mean of errors
% s = standard deviation of errors

function qmeas = generate_discharge_measurements(true_discharge, errortype, mn, s, varargin)

% Set defaults
if nargin>1
    plotflag = varargin{1};
    swot = varargin{2};
else
    plotflag = 0;
    swot = 0;
end

Cv = (s)^2; % variance

switch errortype
    case 'norm'
        disp('Using additive Gaussian errors')
    case 'lnorm'
        disp('Using multiplicative lognormal errors')
end

[nt, m] = size(true_discharge);

% Generate synthetic measurements

mu1 = log((mn^2)/sqrt(Cv+mn^2));
sigma1 = sqrt(log(Cv/(mn^2)+1));
v = lognrnd(mu1, sigma1, 1e5, 1);

qmeas = true_discharge.*lognrnd(mu1, sigma1, nt, m);

% What if we have sparse temporal measurements?

if swot
    true_discharge_SWOT_sampling = NaN(nt,m);
    true_discharge_SWOT_sampling(1:10:end,:) = true_discharge(1:10:end,:);
    qmeas = true_discharge_SWOT_sampling.*lognrnd(mu1, sigma1, nt, m);
end

if plotflag
    figure
    imagescnan(qmeas)
    xlabel('gage')
    ylabel('time')
    title('Discharge measurements')
    colorbar

    figure
    plot(true_discharge)
    hold on
    plot(qmeas, 'r.')
    legend('True discharge','Synthetic measurements')
    xlabel('Time')
    ylabel('Q (mm/day)')
end

return