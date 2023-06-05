% Uses linear routing model to calculate discharge at gages
%
% Uses "dumb" method (not *confusing* like PW13 method)
% 
% INPUTS
% runoff = runoff matrix (nt by n)


function discharge = state_model_dumb(runoff, HH)

disp('Running forward routing model')

% Make robust to NaNs
if any(isnan(runoff(:)))
    % find cells with NaN values
%     sum(isnan(runoff))'
    nancells = find(sum(isnan(runoff))');
    nnan = length(nancells);
    runoff(:, nancells) = [];
    HH(:,nancells,:) = [];
    disp(['Found (and removed) ' num2str(nnan) ' cells with NaN values'])
end

nt = size(runoff,1);
[m,n, kp1] = size(HH);
k = kp1 - 1;

y = NaN(nt, m);
% y = zeros(nt, m);
% i = 0;

x = runoff';

% for tt = ((k+1)+1):nt
for tt = ((k)+1):nt
    %     x = runoff';
    sum1 = 0;
    for ii = 1:(k+1)
        sum1 = sum1 + HH(:,:,ii)*x(:,tt-(ii-1));
    end
    y(tt,:) = sum1;
    
end
discharge = y;

% first (k+1) time steps are spin-up. Can ignore them for now.
quiet = 0;
if ~quiet
    disp('Calculated discharge')
    disp(['First ' num2str(k) ' time steps are spin-up'])
end

return

%% Scratch

%     y(tt, :) = HH(:,:,1)*x(:,tt) + HH(:,:,2)*x(:,tt-1) + HH(:,:,3)*x(:,tt-2);
%     y(tt,:) = HH(:,:,1)*x(:,1) + HH(:,:,2)*x(:,2) + HH(:,:,3)*x(:,3);
%     y(tt,:) = HH(:,:,1)*x(:,1) + HH(:,:,2)*x(:,2);
