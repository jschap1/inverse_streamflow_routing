% Inverse routing
% 
% 6/24/2022
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
%
% 11/12/2023 Attempted to add discharge error tracking (JRS)
% (this feature doesn't work entirely. It produces a jagged pattern
% that repeats every s-k timesteps. Most likely because we aren't
% accounting for Lx_{t-k} cross-correlation terms in the var(y_t)
% calculation. Can look into this in the future.

function [runoff, post_stddev, prior_stddev, K] = ...
    ISR_Y20(runoff_init, HH, gage, s, basin, optionsfile, varargin)

cKF = 0; % flag for constrained KF
pauseloc = 240; % i2 to pause at, if desired

if length(unique(runoff_init(:))) == 1
    % if runoff initial guess is uniform and constant, then K is
    % time-invariant
    uniform_prior = 1;
else
    uniform_prior = 0;
end
% basin.true_runoff = basin.true_runoff';

% flag for using correlation or covariance
usecorr = 0;
if usecorr
    disp('Using correlation instead of covariance')
end
 
% read options
opt = parse_options_file(optionsfile);
opt.quiet = 1;
opt.endwindow = 1e5;

if nargin>6
    [L,T,rho_thres, alpha] = varargin{:};
    opt.L = L;
    opt.T = T;
    opt.rho_thres = rho_thres;
    opt.alpha = alpha;
end

% Compute key parameters
nt = size(gage,1);
[m,n,kp1] = size(HH); % ngages, ncells
k = kp1 - 1; % maximum travel time

% Build augmented state model
[H, L] = calc_augmented_model(HH, s+1);

% Make H and L sparse
H = sparse(H);
L = sparse(L);

% If necessary
H=single(full(H));
L=single(full(L));

% load SC_001.mat
% load TC_001.mat
% rho = SC.*TC;

% Calculate number of smoothing windows
nwin = floor((nt-(k+1)*2)/((s+1)-(k+1))) + 1;

% Handle the steps that aren't fully updated at the end of the ISR
i_last = nwin*(s+1) - (nwin-1)*(k+1);
if ~opt.quiet
    disp(['Last ' num2str(nt - i_last + 1) ' steps will not be fully updated'])
end

runoff = runoff_init;
post_stddev = zeros(nt,n); % standard deviation of runoff estimates
prior_stddev = zeros(nt,n); % standard deviation of runoff prior

% (for future incorporation)
% post_stddevQ = zeros(nt,m); % standard deviation of posterior Q estimates
% prior_stddevQ = zeros(nt,m); % standard deviation of prior Q estimates

if opt.quiet
    1;
else
    tic
end

% keep track of runoff errors
% prior_re = zeros(s+1, n, nwin);
% post_re = zeros(s+1, n, nwin);

always_same_K = 0; % special case for uniform initial runoff
if always_same_K
    load('/Volumes/HD3/ISR/02_Basins/Ohio/PW13_consistency/49gages/assume_L100_T0/kalman_gain_L100_T0.mat')
%     load('./kalman_gain_possibly_illcond.mat')
end
    
% loop over smoothing windows
w = 1; % window number
i2 = 0; % initialize while loop
% for w=1:nwin
while i2<=nt % until out of range
    
    if w==1
        i1 = k+1+1; % start from the first possible window
        i2 = i1 + s; % where x_{t-k} is defined
    end   
    
    % Option to cut off the ISR early
    if w>= opt.endwindow
        break
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
    
    % Get runoff for this window
    x = flipud(runoff_init(i1:i2, :));
    x = reshape(x', n*(s+1), 1); 
    
    % Get x_{t-k} for this window
    if i1-k<=0
        continue
    else
        xtmk = flipud(runoff((i1-k):(i2-k), :));
    end
    xtmk = reshape(xtmk', n*(s+1), 1);
    
    % Get discharge for this window
    y = flipud(gage(i1:i2, :));
    y = reshape(y', m*(s+1), 1);
    
    % Find any missing observations (day without obs)
    missing_rows = isnan(y); 

    % Check for bad values of T and L
    if isinf(opt.T)
        % when T = Inf, P is not invertible
%         opt.T = 9999;
    end
    
    % Decide whether to use true runoff (if available) or prior runoff
    % (realistic scenario) to estimate the Cx matrix
    if opt.use_true_runoff
        xtrue = flipud(basin.true_runoff(i1:i2, :)); % get flipud right...
        xtrue = reshape(xtrue', n*(s+1), 1);         
        x4cov = xtrue;
    else
%         x4cov = ones(size(x)); % 
        x4cov = x;
    end   
    
    % Calculate runoff error covariance
    if ~uniform_prior || ~exist('Cx')
    if ~always_same_K
           
        % can't use this if using domain localization bc domain size
        % changes from subbasin to subbasin
        SC_TC_exist = exist(opt.SC_savename, 'file') > 0 && exist(opt.TC_savename, 'file') > 0;
        if SC_TC_exist % & ~opt.use_domain_loc
            % re-use TC, SC

            %disp(['Calculating prior covariance for window ' num2str(w) ' of ' num2str(nwin)])
            %tic
            [Cx, rho_x] = calc_yang_covariance2(basin, opt, n, s, x4cov, opt.rho_savename);
            %toc

    %         figure
    %         imagesc(Cx(1:5*n,1:5*n))
    %         colorbar

        else
            % calculate TC, SC
            [Cx, rho_x] = calc_yang_covariance2(basin, opt, n, s, x4cov);
%             figure,imagesc(Cx), colorbar, title('Cyy Y20'), set(gca, 'fontsize', 18)
        end       
    end
    
    end
    
    % Calculate runoff error correlation (if desired)
%     SC = load(opt.SC_savename);
%     SC = SC.SC;
%     TC = load(opt.TC_savename);
%     TC = TC.TC;
%     corrmat = SC.*TC;
        
    % Update runoff estimate
    predicted_measurement = H*x + L*xtmk;
    innovation = y - predicted_measurement;
    
    % Construct discharge error covariance matrix
    R = speye(m*(s+1));
    diagR = (opt.sigma^2)*y.^2;
    diagR(missing_rows) = 0; % missing rows could introduce NaNs into R
    R(1:(m*(s+1)+1):end) = diagR;    
    R = single(full(R));
    
  pauseloc = 171;
  if i2>=pauseloc
        1;
  end
    
    % Calculate Kalman gain (using corrmat for exp. purp.)
    if ~always_same_K % in the case where initial runoff is uniform, K will always be the same
    if usecorr
        % use correlation only (ignore covariance)
        Cx = opt.alpha^2*rho_x;
    end
        % use covariance (scaling correlation by runoff initial guess)
%         tic
        [K,H1] = kalman_gain_y20(H, Cx, R, missing_rows, w); % 8 minutes for swot ohio...
%         toc
%     figure, imagesc(corrmat), colorbar, title('\rho')
%     figure, imagesc(K), colorbar, title('K')
    end
    
%     save('./kalman_gain_L100_T0.mat','K','Cx', '-v7.3')

%save('~/Downloads/Y20_true.mat', 'Cx','H','R')    

    % need to get rid of NaN values in innovation to avoid NaNs in x_update
    innovation(missing_rows) = []; % K entries are zero for these rows, so no update
    
    x_update = x + K*innovation;
    Cx_update = (eye(n*(s+1)) - K*H1)*Cx;
    v1 = diag(Cx_update); % posterior runoff variance
    s1 = sqrt(v1); % posterior runoff stddev
    
%     (for future incorporation)
%     posterior discharge variance
%     ptdv = H*Cx_update*H'; 
%     ptdv1 = diag(ptdv); 
%     ptds1 = sqrt(ptdv1);
% 
%     prior discharge variance
%     prdv = H*Cx*H';
%     prds1 = sqrt(diag(prdv));
    
    % constrained KF (optional)
    if cKF
        x_update(x_update<0) = 0;
    end
    
%     if w == 5
%         make_ISR_plot(x, x_update, basin, s, i1, i2)
%     end
    
    % Calculate prior runoff error
%     prior_re(:,:,w) = basin.true_runoff(i1:i2,:) - runoff(i1:i2, :);
    
    % -----Put x back into the runoff array-----
    
    x_update_matrix = flipud(reshape(x_update,n, s+1)');
    stddev_update_matrix_post = flipud(reshape(s1, n, s+1)');
    
    v1 = diag(Cx); % prior variance
    s1 = sqrt(v1); % prior stddev
    stddev_update_matrix_prior = flipud(reshape(s1, n, s+1)');
    
    % (for future incorporation)
%     stddev_update_matrix_postQ = flipud(reshape(ptds1, m, s+1)');
%     stddev_update_matrix_priorQ = flipud(reshape(prds1, m, s+1)');
    
    runoff(i1:i2, :) = x_update_matrix;
    post_stddev(i1:i2, :) = stddev_update_matrix_post;
    prior_stddev(i1:i2, :) = stddev_update_matrix_prior;
    
    % (for future incorporation)
%     post_stddevQ(i1:i2, :) = stddev_update_matrix_postQ;
%     prior_stddevQ(i1:i2, :) = stddev_update_matrix_priorQ;
    
    opt.plot=0;
    if opt.plot
    %%
        % mean posterior discharge
        postypred = state_model_dumb(runoff, HH);

        % mean predicted measurement
        ypred = unfliparoo(predicted_measurement, m, s);

        figure
        gg=1; % gage number
        plot(i1:i2, ypred(:,gg),'blue-*','linewidth',2)
        hold on
        plot(i1:i2, postypred(i1:i2,gg),'red','linewidth',2)
        plot(i1:i2, gage(i1:i2, gg),'.k', 'markersize', 20)
        legend('prior','posterior','measurement')
%         xlim([5,15])
%         ylim([100,250])
        xlabel('timestep')
        ylabel('discharge (mmd)')

        %%
    end

    
    
    
%     
%     x_matrix = unfliparoo(x, n, s);
%     lw = 2;
%     cellnum = 130;
%     figure
%     plot(x_matrix(:,cellnum), 'linewidth', lw)
%     hold on
%     plot(x_update_matrix(:,cellnum), 'linewidth', lw)
%     plot(basin.true_runoff(i1:i2,cellnum), 'linewidth', lw)
%     legend('prior','posterior','truth')
%     
%     % how often is the prior too low?
%     
%     figure
%     plot(mean(x_matrix,2), 'linewidth', lw, 'color', 'blue')
%     hold on
%     plot(mean(x_update_matrix,2), 'linewidth', lw, 'color', 'red')
%     plot(mean(basin.true_runoff(i1:i2,:),2), 'linewidth', lw, 'color', 'black')
%     legend('prior','posterior','truth')

    % Calculate posterior runoff error 
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

return