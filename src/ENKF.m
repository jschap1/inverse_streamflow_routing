% EnKF update function
%
% Modified from Steve's MODWET toolbox
% 
% varargin: opt, rho_yz, rho_zz, rho_yy, H

function [Y_UPDATE]=ENKF(Y_ENS,Z_ENS,ZMEAS,Cv, missing_rows, opt, varargin) 
    % Apply the EnKF to a generic ensemble of states and predicted measurement
    
    if all(missing_rows==1)
        Y_UPDATE = Y_ENS;
        return
    end
    
    if nargin==7
        H = varargin{1};
    elseif nargin==8
        H = varargin{1};
        ind = varargin{2};
    end
    
    % Grab size of ensemble
    Nreps=size(Y_ENS,2);
    % Pre-allocate update 
    Y_UPDATE=NaN(size(Y_ENS));
    
    if opt.multerrors
        CV=diag(Cv*ones(size(Z_ENS,1),1));
    else
        % additive errors
        CV=diag(Cv*ones(size(Z_ENS,1),1));
        m=size(CV,1);
        CV(1:m+1:end) = Cv*ZMEAS.^2;
    end   

    % Compute sample mean vectors
    ymean=mean(Y_ENS,2)*ones(1,Nreps);
    zmean=mean(Z_ENS,2)*ones(1,Nreps);
    zmean(missing_rows,:) = [];
    Z_ENS(missing_rows,:) = [];
    ZMEAS(missing_rows,:) = []; % new
    H(missing_rows,:) = [];
    
    if opt.multerrors
        % using derived form for Cyz and Czz
%         x_ens = exp(Y_ENS);
%         xmean=mean(x_ens,2)*ones(1,Nreps);
%         mu1 = opt.mu1*ones(size(ymean,1),1);
%         mm = mu1*mu1';
%         Cyy=((Y_ENS-ymean)*(Y_ENS-ymean)')/(Nreps-1);
%         Cyz=(Cyy)*H'; % includes bias correction
%         Czz=H*(Cyy)*H';
        Cyz=((Y_ENS-ymean)*(Z_ENS-zmean)')/(Nreps-1);
        Czz=((Z_ENS-zmean)*(Z_ENS-zmean)')/(Nreps-1);
    else
        % additive errors
        Cyz=((Y_ENS-ymean)*(Z_ENS-zmean)')/(Nreps-1);
        Czz=((Z_ENS-zmean)*(Z_ENS-zmean)')/(Nreps-1);
    end
    
    CV(missing_rows,:) = [];
    CV(:,missing_rows) = [];

    % Calculate Kalman gain
    try
        kalm = Cyz/(Czz+CV);
%         kalm = Cyz*pinv(Czz+CV);
    catch
        warning('predicted_meas are nan')
        kalm = zeros(size(Cyz)); 
    end

    % Measurement error realizations needed for update
    v=mvnrnd(zeros(1,size(Czz,1)),CV,Nreps)';
    
    % Additive update
    if sum(missing_rows)==0
        
        ZMEAS(isnan(ZMEAS)) = 0; % so any NaNs don't get in the way...
        for k=1:Nreps
            Y_UPDATE(:,k)=Y_ENS(:,k)+kalm*(ZMEAS+v(:,k)-Z_ENS(:,k));
%             qpost = H*Y_UPDATE(:,k); % this should equal ZMEAS if no measerr
        end
        
    else 
        
        for k=1:Nreps
            Y_UPDATE(:,k)=Y_ENS(:,k)+kalm*(ZMEAS+v(:,k)-Z_ENS(:,k));
        end
        
    end

    opt.plot = 0;
    if opt.plot
        Cyy=((Y_ENS-ymean)*(Y_ENS-ymean)')/(Nreps-1);
        figure,imagesc(Cyy), title('Cyy, Ens'), colorbar, set(gca, 'fontsize', 18)
        figure,imagesc(Cyz), title('Cyz, Ens'), colorbar, set(gca, 'fontsize', 18)
        figure,imagesc(Czz), title('Czz, Ens'), colorbar, set(gca, 'fontsize', 18)
        figure,imagesc(kalm), title('kalm, Ens'), colorbar, set(gca, 'fontsize', 18)
    end

return