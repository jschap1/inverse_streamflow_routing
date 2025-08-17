% Inverse routing
% 
% 12/23/2022
% EnKF version of ISR
%
% INPUTS
% runoff_prior = prior runoff ensemble (n, nt, M)
% HH = measurement model operator
% discharge = gage measurements
% s = length of smoothing window
%
% OUTPUTS
% runoff = posterior runoff ensemble
% w1 = start times for each window
% w2 = end times for each window

function [runoff, w1, w2] = ISR_domain(runoff_prior, HH_sub, gage, s_sub, basin, Cv, opt)

runoff_prior = permute(runoff_prior, [2,1,3]); % put in format (nt, n, M)
small_number = 1e-4; % to avoid taking log of zero (in case there are zeros in the initial guess)

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
[H_sub, L_sub] = calc_augmented_model(HH_sub, s_sub+1);
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

runoff_prior_sub = runoff_prior(:,opt.ind_current_basin,:); % domISR
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
    
    x1_sub = flipud(runoff_prior_sub(i1:i2,:,:));
    x_sub = reshape(permute(x1_sub,[2,1,3]), n_sub*(s_sub+1), M);
    
    % when x is calculated, it should be done using posterior estimates
    % for upstream cells
    x1 = flipud(runoff_prior(i3:i4,:,:));
    x = reshape(permute(x1,[2,1,3]), n*(s+1), M);
    
    % Get x_{t-k} for this window
    if i3-k<=0
        continue
    else
        xtmk = flipud(runoff((i3-k):(i4-k), :,:));
        xtmk_sub = flipud(runoff_sub((i1-k_sub):(i2-k_sub), :,:));
    end
    xtmk = reshape(permute(xtmk,[2,1,3]), n*(s+1), M);
    xtmk_sub = reshape(permute(xtmk_sub,[2,1,3]), n_sub*(s_sub+1), M);

    % Get discharge for this window
    y = flipud(gage(i1:i2, :));
    y = reshape(y', m_sub*(s_sub+1), 1);
    
    % Get discharge for this window, minus any upstream discharge 
    % (the upstream discharge has already been updated)
       
    % Generate the ensemble of discharge errors
    % Assume normally-distributed and iid for now
%     mu_y = zeros(m*(s+1),1);
%     sigma_y = Cv;
%     y_errs = mvnrnd(mu_y, sigma_y, M)';
    

    % to pass H into ENKF (if desired)
    [Htmp, Ltmp] = calc_augmented_model(HH_sub, s_sub+1);
    
    % Predicted measurement for all cells upstream of gauge
    predicted_measurement = H_sub*x + L_sub*xtmk;
    
    % Predicted measurement for cells in subdomain
    predicted_measurement_sub = Htmp*x_sub + Ltmp*xtmk_sub;
    
    % Adjust gauge observation by subtracting contribution (posterior) from cells in
    % upstream subbasins
    opt.subtract_us_runoff = 0;
    if opt.subtract_us_runoff
        tmp1 = ismember((1:n)', opt.ind_current_basin); % get cells that are in upstream basin only
        tmp2 = repmat(tmp1, s+1, 1);
        x_us = x;
        x_us(tmp2,:) = 0;
        xtmk_us = xtmk;
        xtmk_us(tmp2,:) = 0;
        ypred_ups = H_sub*x_us + L_sub*xtmk_us; % predicted meas using cells in upstream subdomains
        y_sub = y - mean(ypred_ups,2);
        if nanmin(y_sub)<0
            warning('negative ysub')
            warning('setting ysub<0 to NaN (skipping update)')
            y_sub(y_sub<0) = NaN;
            % using upstream cells is problematic bc postQ can be larger than
            % trueQ, meaning upstream Q can be larger than downstream Q
        end
    end
            
    missing_rows = isnan(y);
    
    % Adjust gauge observation by subtracting contribution (posterior) from gages in
    % upstream subbasins
    opt.yups=0;
    if opt.yups
        yups = zeros(size(y));
        if min(missing_rows)==0
            nupgages = length(opt.current_upstream_gages);
            if nupgages>0 % if there are any upstream gauges
                for gind=1:nupgages
                    usgagei = opt.current_upstream_gages(gind);
                    lag = opt.us_gages_lag(usgagei);
                    yups1 = opt.gages(i1-lag:i2-lag,usgagei); % streamflow at upstream gauges
                    yups1 = flipud(yups1);
                    yups1 = reshape(yups1', m_sub*(s_sub+1), 1);
                    yups = yups + yups1;
    %                 y_sub = y - yups;
                end
            end
        end
        y_sub = y - yups;
        if min(y_sub)<0
            warning('negative ysub, setting to small number')
            % could also set to NaN to skip
            y_sub(y_sub<0) = small_number;
        end
    end
    
    % Use routing model, but assign discharge to cell
    %ypred_ups = H_sub*x_us + L_sub*xtmk_us;
    
    % Find any missing observations (day without obs)
    missing_rows = isnan(y);
    if min(missing_rows)==0
        1;
    end
    
    if opt.plot
        
        % Compare prior mean runoff to truth
        cellnum = 6;
        mean_prior_runoff_at_grid_cell = flipud(squeeze(mean(x1_sub(:,cellnum,:),3)));
        figure
        plot(i1:i2, mean_prior_runoff_at_grid_cell)    
        hold on
        plot(i1:i2, basin.true_runoff(cellnum, i1:i2), 'k')
        legend('Prior','Truth')
        xlabel('Timestep')
        ylabel('Runoff (mmd)')
        title(['Runoff at cell ' num2str(cellnum)'])
        
        % Compare predicted measurements to measurements (with errors)
        figure
        plot(mean(predicted_measurement,2), 'b-')
        hold on
        plot(y,'r-*')
        legend('Pred. meas','Meas')
        
    end
    
    if i1>40 % 385, 19, 30
        30;
    end
    
    % Use Steve's ENKF function
%     [x_update] = ENKF(x,predicted_measurement,y,Cv, H);   
    
    % Do EnKF update of the log state and observations
    mean1 = 1; % mean of the errors
    mu1 = log((mean1^2)/sqrt(Cv^2+mean1^2)); % mean of the log errors
    sigma1 = sqrt(log(Cv^2/(mean1^2)+1));   
    % ^technically, log(v) ~ N(-0.01, 0.1492) if v ~ LN(1, 0.15), but we'll say 0.15 is close enough    

    % adjust to ensure the measurements are unbiased? (is this necessary?) (no, it is not)
%     tmp = ENKF(log(x), log(predicted_measurement), log(y)-mu1, sigma1^2, opt, rho_yz, rho_zz, missing_rows);

    % Calculate jacobian (experimental)
%     H1 = calculate_jacobian(log(mean(x,2)), log(mean(xtmk,2)), H, L, s, k, 'log2');

    % no log
%     tmp = ENKF(x, predicted_measurement, y, Cv, opt, rho_yz, rho_zz, rho_yy, missing_rows, H);
%     x_update = tmp;
    
%     % eigenvalue decomposition of HH
%     [V,D] = eig(full(H'*H));
%     % V = matrix of eigenvectors
%     % D = diagonal matrix of eigenvalues
%     for j=1:size(D,1)
%         for k=1:j
%             g(j) = D(k,k); % cumulative energy content for each eigvect
%         end
%     end
%     figure,bar(flipud(diag(D)))
%     
%     figure
%     bar(cumsum(flipud(diag(D)))./sum(D(:)))
    
% something is happening to case nan runoff in window 218, subbasin 15
if w==217
    1;
end

    % log
%     small_number = 0;
    if all(missing_rows==1) % debugging
        1;
    elseif sum(missing_rows>0)
        1;
    elseif sum(missing_rows==0)
        1;
    end
        
    if opt.gauss
        % Assuming additive Gaussian error model
        %x_update = ENKF(x, predicted_measurement, y, Cv, missing_rows, Htmp, opt.ind_current_basin);
%         x_update = ENKF(x_sub, predicted_measurement, y, Cv, missing_rows, opt); % old version 
        x_update = ENKF(x_sub, predicted_measurement, y, Cv, missing_rows, opt, Htmp);
        %[Y_UPDATE]=ENKF(Y_ENS,Z_ENS,ZMEAS,Cv, missing_rows, opt, varargin) 
        %x_update = ENKF(x_sub, predicted_measurement_sub, y_sub, Cv, missing_rows, Htmp); % new version
    else
        % Assuming multiplicative lognormal error model
        small_number = 1e-4; % to avoid taking log of zero (in case there are zeros in the initial guess)
        % tmp = ENKF(log(x_sub+small_number), log(predicted_measurement_sub), log(y_sub), Cv, missing_rows, Htmp);
        tmp = ENKF(log(x_sub+small_number), log(predicted_measurement), log(y), Cv, missing_rows, Htmp);
        x_update = exp(tmp) - small_number;
        if min(x_update(:))<0
            warning('negative runoff')
            % x_update(x_update<0) = 0; % in case exp(tmp) = 0
        end
    end
    
    
%     tmp = ENKF(log(x_sub), log(predicted_measurement), log(y), Cv, missing_rows);

    % Ensemble square root filter implementation
    %
    % 8/28/2023 JRS
    % Eliminates the need to simulate observation errors/may improve the update as a result
    
    

    
%     x_update(x_update<0) = 0; % in case exp(tmp) = 0
    
%     postQ = state_model_dumb(x_update(end-71:end,:)', HH_sub);
    
    
%     max(abs(x_update(:) - x_sub(:)))
    
    if max(x_update(:))>1e6 % also chceck for nan/imag
        warning(['Window is ' num2str(w)])
        warning('x is too big')
    end
    
%     xmean = mean(x,2)*ones(1,M);
%     ymean = mean(y + y_errs,2)*ones(1,M);
%     Cxy = ((x-xmean)*((y + y_errs)-ymean)')/(M-1);
%     Cyy =(((y + y_errs)-ymean)*((y + y_errs)-ymean)')/(M-1);
%     K1 = Cxy/[Cyy+Cv];
    
%     % Estimate Kalman gain
%     C = cov(x');
%     figure, imagesc(C), title('Estimated covariance'), colorbar
%     K = C*H'*inv(H*C*H' + Cv);
%     figure, imagesc(K), title('Estimated K'), colorbar
    
%     x_update = x + K*innovation;

    % Put x back into the runoff array
    % (need to get indexing right!!)
        
    x_update_matrix = unfliparoom(x_update, n_sub, s_sub); 
    runoff(i1:i2,opt.ind_current_basin,:) = x_update_matrix;
%     aa = reshape(x_update, s+1, n, M);
%     x_update_matrix = flipud(a);
%     ro = runoff(:,opt.ind_current_basin,:); 
%     postypred = state_model_dumb(mean(ro,3), HH_sub);
%     if min(postypred(:))<0
%         1
%     end
if opt.plot
%%
    % mean posterior discharge
    ro = runoff(:,opt.ind_current_basin,:); 
%     postypred = state_model_dumb(mean(runoff,3), HH);
    postypred = state_model_dumb(mean(ro,3), HH_sub);
    
    % mean predicted measurement
    ypred = mean(unfliparoom(predicted_measurement_sub, m_sub, s_sub),3);
    
    figure
    gg=1; % gage number
    plot(i1:i2, ypred(:,gg),'blue-*','linewidth',2)
    hold on
    plot(i1:i2, postypred(i1:i2,gg),'red','linewidth',2)
    plot(i1:i2, flipud(y_sub),'.k', 'markersize', 20)
%     plot(i1:i2, gage(i1:i2, gg),'.k', 'markersize', 20)
    legend('prior','posterior','measurement')
    xlim([10,50])
    ylim([0,100])
    xlabel('timestep')
    ylabel('discharge (mmd)')
    
end

%     cell6_updated_mean_runoff = mean(x_update_matrix(:,6,:),3);
       
    if opt.plot
    if i1==385
        ms = 16;
        figure
%         plot(basin.t(i1:i2), gage(i1:i2,5)) 
        plot(basin.t(i1:i2), basin.true_discharge(i1:i2,5)) 
        hold on
        plot(basin.t(i1:i2), gage(i1:i2,5), '.', 'markersize', ms)
        legend('Truth','Measurements')
        title('Discharge at Gage 5')
    end
    end
    
    % Calculate posterior runoff error
%     Cpost = (eye(n*(s+1)) - K*H)*C;   
%     runoff_err = unfliparoo(diag(Cpost),n,s)'; % get the variance
%     posterior_runoff_variance(:, i1:i2) = runoff_err; % put in array
%     posterior_runoff_variance(:, i1:i2) = zeros(n, s+1); % filler

    % Convince self that my method works the same as Steve's
    % (using C = cov(x);
    
%     post_re(:,:,w) = basin.true_runoff(i1:i2,:) - runoff(i1:i2, :);
    
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

runoff = permute(runoff, [2,1,3]);

% Check if posterior discharge matches observed discharge
% opt.ind_current_basin

% prior.mean_runoff = mean(runoff_prior,3);
% posterior.mean_runoff = mean(runoff,3);
% 
% prior.meanQ = state_model_dumb(prior.mean_runoff(:,opt.ind_current_basin), HH_sub);
% posterior.meanQ = state_model_dumb(posterior.mean_runoff(opt.ind_current_basin,:)', HH_sub);
% 
% figure
% plot(prior.meanQ, 'b', 'linewidth', 2)
% hold on
% plot(posterior.meanQ, 'r', 'linewidth', 2)
% plot(gage, 'k.', 'markersize', 20)
% legend('Prior Q','Post Q', 'Observations')


return