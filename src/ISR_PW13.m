% Inverse routing
% 
% 6/24/2022
% Simplest case: PW13
%
% INPUTS
% initial_runoff
% HH = measurement model operator
% discharge = gage measurements
% s = length of smoothing window
% a = alpha^2
% R = measurement error covariance matrix
%
% OUTPUTS
% runoff = updated runoff estimate
% [w1,w2] beginning and end of each window

function [runoff, K, H, w1, w2] = ISR_PW13(runoff_init, HH, gage, s, P_option, alpha1, meas_err, varargin)

opt.endwindow = 1e10;

if nargin==8
    basin = varargin{1};
    basinflag = 1;
else
    basinflag = 0;
end

opt.timer = 1; tic
opt.progress = 1;
opt.plotP = 0;

% Compute key parameters
nt = size(gage,1);
[m,n,kp1] = size(HH); % ngages, ncells
k = kp1 - 1; % maximum travel time

% Build augmented state model
[H, L] = calc_augmented_model(HH, s+1);

% Calculate number of smoothing windows
nwin = floor((nt-(k+1)*2)/((s+1)-(k+1))) + 1;

% Handle the steps that aren't fully updated at the end of the ISR
% i_last = nwin*(s+1) - (nwin-1)*(k+1);
% disp(['Last ' num2str(k+1) ' steps will not be fully updated'])

runoff = runoff_init;

% For efficiency
H = sparse(H);
L = sparse(L);
% R = sparse(R);

no_missing_data = 1; % flag for the no_missing_data case (daily obs)
disp('Check no_missing_data flag')
proportional = 0;
constK = 0;
switch P_option
    case 'const_diag'
        % Under this assumption, K does not change from window to window
        P = alpha1^2*speye(n*(s+1));
        if meas_err == 0
            constK = 1;
            % K is only truly constant if there is no measurement error
            if no_missing_data == 1
                K = H'/(H*H');
            end
            % Also, if there are missing rows, K will have to be
            % recalculated
        end   
    case 'proportional'
        P = alpha1^2*speye(n*(s+1));
        proportional = 1;
    case 'filled'
        P = a*ones(n*(s+1),n*(s+1)); 
    case 'negative'
        P = a*eye(n*(s+1));
        % ensures the runoff error in cells 1 and 2 are always negatively
        % correlated in space (with correlation coef = 1)
        b = 0.5;
        P(2:7:end) = -1*b;
        P(7:7:end) = -1*b;
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
    
    w1(w) = i1;
    w2(w) = i2;
    
    if w>= opt.endwindow
        break
    end    
    
    % we add an extra +1 to i1 in case discharge is 0 for the first k+1
    % time steps, as it is for the 2cell synthetic experiments, when the
    % routing model is used to generate true discharge from true runoff -
    % the routing model has k+1 timesteps of spinup
    
%     i1 = 1 + (w-1)*(s-k);
%     i2 = w*(s+1) - (w-1)*(k+1); 
    
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
    
    % Calculate runoff error covariance
    if proportional
%         alpha1 = a;
        P(1:(n*(s+1)+1):end) = alpha1^2*x.^2;
        if basinflag
            xtrue = fliparoo(basin.true_runoff(:,i1:i2)', s);
            P(1:(n*(s+1)+1):end) = alpha1^2*xtrue.^2;
        end
    end
    
    % Update runoff estimate
    predicted_measurement = H*x + L*xtmk;
    innovation = y - predicted_measurement;
    
%     figure
%     plot(y)
%     hold on 
%     plot(predicted_measurement)
%     legend('true','pred')

    % Construct measurement error matrix
    R = speye(m*(s+1));
%     diagR = ones(m*(s+1),1)*meas_err^2;
    diagR = meas_err*y.^2;
    diagR(missing_rows) = 0; % missing rows could introduce NaNs into R
    R(1:(m*(s+1)+1):end) = diagR;
    
    if strcmp(P_option, 'proportional') || (strcmp(P_option, 'const_diag') && ~constK || sum(missing_rows)>0)
        K = kalman_gain_pw13(H, P, R, missing_rows, w);
%         figure,imagesc(P(1:n,1:n)),colorbar, title('P')
%         figure,imagesc(R(1:m,1:m)),colorbar, title('R')
%         figure,imagesc(K(1:n,1:m)),colorbar, title('K')
        % otherwise K could be time-invariant
    end
    
    % Alternative Kalman gain formulation (still working on this one)
%     for i=1:n(s+1)
%         for j=1:m(s+1)
%             K(i,j) = x(i)/sum(x(k)^2);
%         end
%     end
    
    if sum(missing_rows)>0
        innovation(missing_rows) = [];
    end
    x_update = x + K*innovation; % by having K as sparse, it means 0*nan = 0
%     P_update = (eye(n*(s+1)) - K*H)*P;

    if any(x_update<0) % any(isnan(x_update))
        sign(mean(innovation));
    end
    
    if opt.plotP
        figure
        imagesc(P)
        title('Runoff error covariance')
        colorbar
        
        figure,imagesc(eye(6)*K), title('Kalman gain'), colorbar, set(gca, 'fontsize', 16), caxis([0,0.5])
        
    end
    
    % Put x back into the runoff array
    x_update_matrix = flipud(reshape(x_update,n, s+1)');
    runoff(i1:i2, :) = x_update_matrix;
    
    % Move to next window
    w = w + 1;
    i1 = i2-k;
    i2 = i1 + s; % ensure k+1 overlapping steps
    
    % Displays time to go
    disp(['Processed window ' num2str(w) ' of ' num2str(nwin)])
    if opt.timer == 1
        time_for_iteration = toc;
        opt.timer = 0;
        if opt.progress
            disp(['Approximately ' num2str(time_for_iteration*nwin) ' seconds to go'])
        end
    end    
    
end % end main loop

% runoff = runoff';

return