% Inverse routing
% 
% 7/21/2023
% Y20 case - spatiotemporal error correlations
%
% INPUTS
% initial_runoff
% HH = measurement model operator
% discharge = gage measurements
% s = length of smoothing window
%
% OUTPUTS
% runoff = updated runoff estimate

function [runoff, w1, w2] = ISR_Y20_domain(runoff_prior, HH_sub, gage, s_sub, basin, opt, L, T, rho_thres)

runoff_prior = runoff_prior; % put in format (nt, n, M)

checknan = sum(isnan(runoff_prior(:)));
if checknan>0
    warning(['Filled ' num2str(checknan) ' missing values in runoff prior'])
    runoff_prior = fillmissing(runoff_prior, 'linear');
    runoff_prior(runoff_prior<0) = 0;
end

% if opt.loc
%     disp('using cov localization')
% end
opt.quiet = 1;
opt.timer = 1;
opt.plot = 0;
opt.progress = 1;
HH = opt.HH_total;

% Compute key parameters
nt = size(gage,1);
[m_sub,n_sub,kp1_sub] = size(HH_sub); % ngages, ncells
k_sub = kp1_sub - 1; % maximum travel time
s = opt.s_total;
[m,n,kp1] = size(opt.HH_total); % ngages, ncells
k = kp1 - 1;

% Build augmented state model
[H_subK, L_subK] = calc_augmented_model(HH_sub, s_sub+1); % used for Kalman gain
[H, L] = calc_augmented_model(HH, s+1);

% we only want H for the current gage and the subdomain's s (from t-s_sub to t)
H_sub = H(1:m*(s_sub+1),:);
H_sub = H_sub(opt.current_gage:m:end,:); 
% need to do something similar for L, but account for the different t-k
L_sub = L(1:m*(s_sub+1),:);
L_sub = L_sub(opt.current_gage:m:end,:); 

% Make H and L sparse (can do for H_sub, too, if desired)
% H = sparse(H);
% L = sparse(L);

% Calculate number of smoothing windows
nwin = floor((nt-(k_sub+1)*2)/((s_sub+1)-(k_sub+1))) + 1;

% Handle the steps that aren't fully updated at the end of the ISR
i_last = nwin*(s_sub+1) - (nwin-1)*(k_sub+1);
disp(['Last ' num2str(nt - i_last + 1) ' steps will not be fully updated'])

M = size(runoff_prior,3); % ensemble size

runoff_prior_sub = runoff_prior(:,opt.ind_current_basin); % domISR
runoff = runoff_prior; % posterior runoff ensemble
runoff_sub = runoff_prior_sub; % initializing for later use (xtmk)

% loop over smoothing windows
start_flag = 1;
w = 0; % window number
i2 = 0; % initialize while loop
% for w=1:nwin
while i2<=nt % until out of range

    w = w + 1;
    
    if start_flag
        i1 = k_sub+1+1; % start from the first possible window
        i2 = i1 + s_sub; % where x_{t-k} is defined
        start_flag = 0;
        
        % also define indices for the full domain
        i4 = i2; % not just the subbasin
        i3 = i4-s;
    else
        % Move to next window
        i1 = i2-k_sub;
        i2 = i1 + s_sub; % ensure k+1 overlapping steps        
        i4 = i2;
        i3 = i4 - s;
    end
    
    % we add an extra +1 to i1 in case discharge is 0 for the first k+1
    % time steps, as it is for the 2cell synthetic experiments, when the
    % routing model is used to generate true discharge from true runoff -
    % the routing model has k+1 timesteps of spinup

% for w=1:nwin

%     i1 = 1 + (w-1)*(s-k);
%     i2 = w*(s+1) - (w-1)*(k+1);

    % Check 
    if i2>nt
        if ~opt.quiet
            disp('Cannot go further, next window is beyond range')
        end
        break
    end
    
    % Check: cannot start until i3>=1
    if i3-k<=0
        continue
    end

    % track windows for output
    w1(w) = i1;
    w2(w) = i2;
    
    % Get runoff for this window
    x1_sub = flipud(runoff_prior_sub(i1:i2,:));
    x_sub = reshape(x1_sub', n_sub*(s_sub+1), 1);
    
    % when x is calculated, it should be done using posterior estimates
    % for upstream cells
    x1 = flipud(runoff_prior(i3:i4,:));
    x = reshape(x1', n*(s+1), 1);
    
    % Get x_{t-k} for this window
    if i3-k<=0
        continue
    else
        xtmk = flipud(runoff((i3-k):(i4-k),:));
        xtmk_sub = flipud(runoff_sub((i1-k_sub):(i2-k_sub),:));
    end
    xtmk = reshape(xtmk', n*(s+1), 1);
    xtmk_sub = reshape(xtmk_sub', n_sub*(s_sub+1), 1);

    % Get discharge for this window
    y = flipud(gage(i1:i2, :));
    y = reshape(y', m_sub*(s_sub+1), 1);
    
    % Find any missing observations (day without obs)
    missing_rows = isnan(y);
   
    % Generate the ensemble of discharge errors
    % Assume normally-distributed and iid for now
%     mu_y = zeros(m*(s+1),1);
%     sigma_y = Cv;
%     y_errs = mvnrnd(mu_y, sigma_y, M)';
    
    % Calculate predicted measurements and innovation
    predicted_measurement = H_sub*x + L_sub*xtmk;
        
    % Calculate covariance
    opt.RAM = 36;
    opt.SC_savename = './tempSC.txt';
    opt.TC_savename = './tempTC.txt';
    opt.rho_savename = './temprho.txt';
    opt.alpha = 1;
    Cx = calc_yang_covariance2(basin, opt, n_sub, s_sub, x_sub);

    innovation = y - predicted_measurement;
    R = sparse(m_sub*(s_sub+1), m_sub*(s_sub+1));
    R(1:(m_sub*(s_sub+1)+1):numel(R)) = opt.sigma^2;
    R = single(full(R));
    
    K = kalman_gain_y20(H_subK, Cx, R, missing_rows, w); % 8 minutes for swot ohio...
    innovation(missing_rows) = 1; % K entries are zero for these rows, so no update
    x_update = x_sub + K*innovation;
    
    % Put x back into the runoff array
    x_update_matrix = flipud(reshape(x_update,n_sub, s_sub+1)');
    runoff(i1:i2, opt.ind_current_basin) = x_update_matrix;
          
%     x_update_matrix = unfliparoom(x_update, n_sub, s_sub); 
%     runoff(i1:i2,opt.ind_current_basin,:) = x_update_matrix;
%     aa = reshape(x_update, s+1, n, M);
%     x_update_matrix = flipud(a);
   
    % Displays time to go
    if ~opt.quiet
        disp(['Processed window ' num2str(w) ' of ' num2str(nwin)])
        if opt.timer == 1
            time_for_iteration = toc;
            opt.timer = 0;
            if opt.progress
                disp(['Approximately ' num2str(time_for_iteration*nwin) ' seconds to go'])
            end
        end    
    end

end % end main loop

return