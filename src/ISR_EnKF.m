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

function [runoff, w1, w2] = ISR_EnKF(runoff_prior, HH, gage, s, basin, Cv, opt)

runoff_prior = permute(runoff_prior, [2,1,3]); % put in format (nt, n, M)

if opt.loc
    disp('using cov localization')
end
opt.quiet = 0;
opt.timer = 1;
opt.plot = 0;
opt.progress = 1;
opt.HH = HH;

% Compute key parameters
nt = size(gage,1);
[m,n,kp1] = size(HH); % ngages, ncells
k = kp1 - 1; % maximum travel time

% Build augmented state model
[H, L] = calc_augmented_model(HH, s+1);

% Make H and L sparse
H = sparse(H);
L = sparse(L);

% Calculate number of smoothing windows
nwin = floor((nt-(k+1)*2)/((s+1)-(k+1))) + 1;

% Handle the steps that aren't fully updated at the end of the ISR
i_last = nwin*(s+1) - (nwin-1)*(k+1);
disp(['Last ' num2str(nt - i_last + 1) ' steps will not be fully updated'])

M = size(runoff_prior,3);
    
% Calculate space-time correlation matrices for covariance localization
% [rho_yy, rho_yz, rho_zz] = calc_localization_matrices(basin, s, opt);
% ^this might change if different gauges are active at diff times?  

% figure,imagesc(rho_zz), title('\rho_{zz}'), colorbar
% figure,imagesc(rho_yz), title('\rho_{yz}'), colorbar
% figure,imagesc(rho_yy), title('\rho_{yy}'), colorbar

% if opt.domISR
%     runoff_prior = runoff_prior(:,opt.ind_current_basin,:); % domISR
% end

runoff = runoff_prior; % posterior runoff ensemble

% loop over smoothing windows
w = 1; % window number
i2 = 0; % initialize while loop
% for w=1:nwin
while i2<=nt % until out of range

    if w==1
        i1 = k+1+1; % start from the first possible window
        i2 = i1 + s; % where x_{t-k} is defined
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

    % track windows for output
    w1(w) = i1;
    w2(w) = i2;
    
    % Get runoff for this window
    x1 = flipud(runoff_prior(i1:i2,:,:));
    x = reshape(permute(x1,[2,1,3]), n*(s+1), M);
    
    % Get x_{t-k} for this window
    if i1-k<=0
        continue
    else
        xtmk = flipud(runoff((i1-k):(i2-k), :,:));
    end
    xtmk = reshape(permute(xtmk,[2,1,3]), n*(s+1), M);

    % Get discharge for this window
    y = flipud(gage(i1:i2, :));
    y = reshape(y', m*(s+1), 1);
    
    % Find any missing observations (day without obs)
    missing_rows = isnan(y);
   
    % Generate the ensemble of discharge errors
    % Assume normally-distributed and iid for now
%     mu_y = zeros(m*(s+1),1);
%     sigma_y = Cv;
%     y_errs = mvnrnd(mu_y, sigma_y, M)';
    
    % Calculate predicted measurements and innovation
    predicted_measurement = H*x + L*xtmk;
%     innovation = y + y_errs - predicted_measurement;
    
    % Pause and see what is going on in the ISR update
            
    if opt.plot
        
        % Compare prior mean runoff to truth
        cellnum = 6;
        mean_prior_runoff_at_grid_cell = flipud(squeeze(mean(x1(:,cellnum,:),3)));
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
        plot(mean(predicted_measurement,2))
        hold on
        plot(y)
        legend('Pred. meas','Meas')
        
    end
    
    if i1==30 % 385, 19, 30
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
    
    % log

%     tmp = ENKF(log(x), log(predicted_measurement), log(y), Cv, missing_rows);
%     tmp = ENKF(log(x), log(predicted_measurement), log(y), Cv, missing_rows, opt, rho_yz, rho_zz, rho_yy, H);

    small_number = 1e-4; % to avoid taking log of zero (in case there are zeros in the initial guess)
    tmp = ENKF(log(x+small_number), log(predicted_measurement), log(y), Cv, missing_rows);
    x_update = exp(tmp) - small_number;
    
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
        
    x_update_matrix = unfliparoom(x_update, n, s); 
    runoff(i1:i2,:,:) = x_update_matrix;
%     aa = reshape(x_update, s+1, n, M);
%     x_update_matrix = flipud(a);
    
if opt.plot
    mxum = mean(x_update_matrix,3);
    figure
    subplot(1,2,1)
    imagesc(mean(runoff(i1:i2,:,:),3)), colorbar, title('prior')
    subplot(1,2,2)
    imagesc(mxum), colorbar, title('update')
end

if opt.plot

    % mean posterior discharge
    postypred = state_model_dumb(mean(runoff,3), HH);
    
    % mean predicted measurement
    ypred = mean(unfliparoom(predicted_measurement, m, s),3);
    
    figure
    gg=1; % gage number
    plot(i1:i2, ypred(:,gg),'blue','linewidth',2)
    hold on
    plot(i1:i2, postypred(i1:i2,gg),'red','linewidth',2)
    plot(i1:i2, gage(i1:i2, gg),'.k', 'markersize', 20)
    legend('prior','posterior','measurement')
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

    % Move to next window
    w = w + 1;
    i1 = i2-k;
    i2 = i1 + s; % ensure k+1 overlapping steps

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

return