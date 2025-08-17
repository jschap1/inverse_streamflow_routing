% EnKF update function
%
% Modified from Steve's MODWET toolbox
% 
% varargin: opt, rho_yz, rho_zz, rho_yy, H

function [Y_UPDATE]=ENKF(Y_ENS,Z_ENS,ZMEAS,Cv, missing_rows, varargin) 
    % Apply the EnKF to a generic ensemble of states and predicted measurement
    
    if all(missing_rows==1)
        Y_UPDATE = Y_ENS;
        return
    end
    
    if nargin==6
        H = varargin{1};
    elseif nargin==7
        H = varargin{1};
        ind = varargin{2};
    end
    
    % Grab size of ensemble
    Nreps=size(Y_ENS,2);
    % Pre-allocate update 
    Y_UPDATE=NaN(size(Y_ENS));
    
    % For multiplicative errors, CV is a percent
    % Define measurement error covariance assuming homoscedastic error with 
    % provided Cv as the diagonal (variance)
%     CV=diag(Cv*ones(size(Z_ENS,1),1));
    
    % For no-log/additive error case, 15% err is var(y) = (0.15*y)^2
    
    CV=diag(Cv*ones(size(Z_ENS,1),1));
    m=size(CV,1);
    CV(1:m+1:end) = Cv*ZMEAS.^2;
    
    % Compute sample mean vectors
    ymean=mean(Y_ENS,2)*ones(1,Nreps);
    zmean=mean(Z_ENS,2)*ones(1,Nreps);
    % Compute sample covariance matrices
%     Cyy=((Y_ENS-ymean)*(Y_ENS-ymean)')/(Nreps-1);
%     Cyz=((Y_ENS-ymean)*(Z_ENS-zmean)')/(Nreps-1);
%     Czz=((Z_ENS-zmean)*(Z_ENS-zmean)')/(Nreps-1);
    
    % What if we account for missing measurements earlier on? yes, good, do
    % this
    zmean(missing_rows,:) = [];
    Z_ENS(missing_rows,:) = [];
    ZMEAS(missing_rows,:) = []; % new
    Cyz=((Y_ENS-ymean)*(Z_ENS-zmean)')/(Nreps-1);
    
    % do stuff
%     A = (Z_ENS-zmean)*(Z_ENS-zmean)';
%     for ii=1:size(A,1)
%         mean(A,1)
%     end
    
    Czz=((Z_ENS-zmean)*(Z_ENS-zmean)')/(Nreps-1);
    CV(missing_rows,:) = [];
    CV(:,missing_rows) = [];
%     kalm = Cyz/(Czz+CV);
%     Czz=eye(size(Czz,1));
%     Y20 = load('~/Downloads/Y20_true.mat', 'Cx','H','R');
%     
%     % gives different results (missing after vs. missing before)
%     Y20.Cyz = Y20.Cx*Y20.H';
%     Y20.Czz = Y20.H*Y20.Cx*Y20.H';
%     Y20.kalm = Y20.Cyz/Y20.Czz;
%     missing_rows(1:7)=0;
%     Y20.kalm(:,missing_rows) = [];

    % gives different results
%     Y20.H1 = Y20.H;
%     Y20.H1(missing_rows,:) = [];
%     Y20.R(missing_rows,:) = [];
%     Y20.R(:,missing_rows) = [];
%     Y20.K = (Y20.Cx*Y20.H1')*pinv(Y20.H1*Y20.Cx*Y20.H1' + Y20.R);
%     figure,imagesc(Y20.K), title('kalm, Y20'), colorbar, set(gca, 'fontsize', 18)
    

% figure,imagesc(Y20.Cyz), title('Cyz, Y20'), colorbar, set(gca, 'fontsize', 18), caxis([0,150])
% figure,imagesc(Y20.Czz), title('Czz, Y20'), colorbar, set(gca, 'fontsize', 18), caxis([0,2.5e4])
% figure,imagesc(Y20.Cx), title('Cyy, Y20'), colorbar, set(gca, 'fontsize', 18), caxis([0,4])
% figure,imagesc(Y20.kalm), title('kalm, Y20'), colorbar, set(gca, 'fontsize', 18)

% figure,imagesc(Cyy), title('Cyy, Ens'), colorbar, set(gca, 'fontsize', 18), caxis([0,4])
% % figure,imagesc(Cyy/4), title('corr(y), Ens'), colorbar

%Czz(Czz<1e-2)=0;
% Cyy(Cyy/4<exp(-2)) = 0; % equivalent localization to Y20
% H1 = H; H1(missing_rows,:) = [];
% Czz = H1*Cyy*H1';
% Cyz = Cyy*H1';

% calculate the correlation matrix
% sqrt(diag(Cyy))

% Cyy(Cyy<0.4)=0;

% Removing rows and columns from Czz and Cyz doesn't help
% Something about the log update is not optimal
% val=21;
% Czz = Czz(1:val,1:val);
% Cyz = Cyz(:,1:val);
% CV = CV(1:val,1:val);
% ZMEAS = ZMEAS(1:val);
% Z_ENS = Z_ENS(1:val,:);

    % Calculate Kalman gain
    try
        kalm = Cyz/(Czz+CV);
%         kalm = Cyz*pinv(Czz+CV);
    catch
        warning('predicted_meas are nan')
        kalm = zeros(size(Cyz)); 
    end

% figure,imagesc(Cyz), title('Cyz, Ens'), colorbar, set(gca, 'fontsize', 18)%, caxis([0,150])
% figure,imagesc(Czz), title('Czz, Ens'), colorbar, set(gca, 'fontsize', 18)%, caxis([0,2.5e4])
% figure,imagesc(inv(Czz)), title('inv Czz, Ens'), colorbar, set(gca, 'fontsize', 18)%, caxis([0,2.5e4])
% figure,imagesc(kalm), title('kalm, Ens'), colorbar, set(gca, 'fontsize', 18), %caxis([-0.10,0.12])


%     if sum(missing_rows)>0 % account for missing rows
%         kalm(:,missing_rows) = [];
%     end            

    % Measurement error realizations needed for update
    v=mvnrnd(zeros(1,size(Czz,1)),CV,Nreps)';
    
    % Additive update
    if sum(missing_rows)==0
        
        ZMEAS(isnan(ZMEAS)) = 0; % so the NaNs don't get in the way...
        Y_CHANGE = Y_UPDATE;
        for k=1:Nreps
            Y_UPDATE(:,k)=Y_ENS(:,k)+kalm*(ZMEAS+v(:,k)-Z_ENS(:,k));
        end
        
    else 
        
        %ZMEAS(missing_rows,:) = []; % new
        %Z_ENS(missing_rows,:) = []; % new
        %v(missing_rows,:) = []; % new
        Y_CHANGE = Y_UPDATE;
        for k=1:Nreps
            Y_UPDATE(:,k)=Y_ENS(:,k)+kalm*(ZMEAS+v(:,k)-Z_ENS(:,k));
        end
        
%         [~,maxi]=max(Y_UPDATE(:))
%         [ii,jj] = ind2sub(size(Y_ENS),maxi)
%         Y_UPDATE(955,754)
%         xupdate = exp(mean(Y_UPDATE,2));
%         y = H*x_update;
        
    end

return