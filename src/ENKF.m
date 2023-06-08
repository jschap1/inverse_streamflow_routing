% EnKF update function
%
% Modified from Steve's MODWET toolbox
% 
% Can do a multiplicative update or an additive update
% Localization
%
% varargin: opt, rho_yz, rho_zz, rho_yy, H

function [Y_UPDATE]=ENKF(Y_ENS,Z_ENS,ZMEAS,Cv, missing_rows, varargin) 
    % Apply the EnKF to a generic ensemble of states and predicted measurement
    
    if nargin>5
        opt = varargin{1};
        rho_yz = varargin{2};
        rho_zz = varargin{3};
        rho_yy = varargin{4};
        H = varargin{5};
        opt.loc = 1;
    else
        opt.loc = 0;
    end
    
    % Grab size of ensemble
    Nreps=size(Y_ENS,2);
    % Pre-allocate update 
    Y_UPDATE=NaN(size(Y_ENS));
    % Define measurement error covariance assuming homoscedastic error with 
    % provided Cv as the diagonal (variance)
    CV=diag(Cv*ones(size(Z_ENS,1),1));
    
    % Account for missing observations (by inflating Cv for missing rows)
    % this may be causing the weird results. See if we can get rid of this
    % and handle missing observations by modifying the Kalman gain matrix
    % directly.
%     if any(missing_rows)
%         lcv = size(CV,1);
%         inflate_vals = Cv*ones(lcv,1);
%         inflate_vals(missing_rows) = 1e9;
%         CV(1:lcv+1:end) = inflate_vals;
%         ZMEAS(missing_rows) = nanmean(ZMEAS); % any value will do
%     end

    % Compute sample mean vectors
    ymean=mean(Y_ENS,2)*ones(1,Nreps);
    zmean=mean(Z_ENS,2)*ones(1,Nreps);
    % Compute sample covariance matrices
    Cyz=((Y_ENS-ymean)*(Z_ENS-zmean)')/(Nreps-1);
    Czz=((Z_ENS-zmean)*(Z_ENS-zmean)')/(Nreps-1);
%     Cyy=((Y_ENS-ymean)*(Y_ENS-ymean)')/(Nreps-1);

%     if opt.loc
%         Cyz_loc = rho_yy.*Cyy*H';
%         Czz_loc = H*(rho_yy.*Cyy)*H';
%         kalm_loc = Cyz_loc/(Czz_loc+CV);
%         if sum(missing_rows)>0 % account for missing rows
%             kalm_loc(:,missing_rows) = 0;
%         end               
%     else
%         Cyz = Cyy*H';
%         Czz = H*Cyy*H';
%         kalm = Cyz*inv(Czz);
%     end
    
%     figure
%     subplot(1,2,1)
%     imagesc(inv(Czz)), colorbar, title('Cyz')
%     subplot(1,2,2)
%     imagesc(inv(Czz_loc)), colorbar, title('Czz')
%     subplot(1,2,2)
%     imagesc(rho_zz), colorbar, title('\rho_{zz}')
    
%     Czz_true = [2,0,0;0,2,0;1,0,1];
%     Czz=[2.01,0.01,0.02;0.05,2.1,0.1;1.1,0.2, 0.98];
%     rho_zz=[1,0.01,0.01;0.01,1,0.01;1,0.01,1];

%     if opt.rho_yy
%         % do localization with rho_yy and compare
%         Cyy=((Y_ENS-ymean)*(Y_ENS-ymean)')/(Nreps-1);
%         Cyy_loc = rho_yy.*Cyy;
%         
%         % These operations are not valid because they are only true for y,z
%         % not for log(y) and log(z) as they are being applied here for
%         % ensemble ISR
%         Cyz = Cyy*H';
%         Czz = H*Cyy*H';        
%         Cyz_loc = Cyy_loc*H';
%         Czz_loc = H*Cyy_loc*H';
        
%         Czz = log(1 + H*Cyy*H'); % need denominator, need to derive something
%         Cyz = log(1+Cyy*H');
%         figure,imagesc(inv(Czz+CV)), colorbar, colormap(bluewhitered)
        
%         kalm = Cyz/(Czz+CV);
%         kalm_loc = Cyz_loc/(Czz_loc+CV);
%         figure
%         subplot(2,4,1), imagesc(Cyy), title('C_{yy}'), colorbar
%         subplot(2,4,2), imagesc(Czz), title('Cz'), colorbar, colormap(bluewhitered), caxis([0,1])
%         subplot(2,4,3), imagesc(inv(Czz + CV)), title('inv Cz'), colorbar, colormap(bluewhitered)
%         subplot(2,4,4), imagesc(kalm), title('K'), colorbar, colormap(bluewhitered)
%         subplot(2,4,5), imagesc(Cyy_loc), title('C_{yy} loc'), colorbar, colormap(bluewhitered)
%         subplot(2,4,6), imagesc(Czz_loc), title('Cz loc'), colorbar, colormap(bluewhitered), caxis([0,1])
%         subplot(2,4,7), imagesc(inv(Czz_loc + CV)), title('inv Cz loc'), colorbar, colormap(bluewhitered)
%         subplot(2,4,8), imagesc(kalm_loc), title('K loc'), colorbar
        
%         Cy = eye(480);
%         Cy(1,1)=10;
%         Cz = H*Cy*H';
%         figure,subplot(1,2,1),imagesc(Cz-Cz10), title('Cz'),subplot(1,2,2),imagesc(Cz10), title('Cz10')
%     end
    
    if opt.loc
        rho_zz_orig = rho_zz;
    end
    if opt.loc & opt.rho_yy==0
        
        % ensures that localization does not remove built-in correlations
        % that come from the routing model - we want to keep those
        yz_mask = H'>0;
%         rho_yz(yz_mask) = 1;
%         rho_yz_SPD = nearestSPD(rho_yz);
        Cyz_loc = rho_yz.*Cyz;
%         Cyz_loc = rho_yz.*Cyz;
        
        zz_mask = H*H'>0;
%         rho_zz(zz_mask) = 1;
%         rho_zz_SPD = nearestSPD(rho_zz);
%         force_symmetric_diagonally_dom(rho_zz)
        Czz_loc = rho_zz.*Czz;
        Czz_loc_orig = rho_zz_orig.*Czz;
        
        % special case
%         Cyz_loc = Cyz;
%         Czz_loc = Czz;
%         Cyz_loc(~yz_mask)=0;
%         Czz_loc(~zz_mask)=0;
        
%         Czz_loc(Czz_loc<1e-1) = 0;
        
%         n=20;
%         m = 5;
%         Cz = Czz(1:m,1:m);
%         Cz_loc = rho_zz(1:m,1:m).*Cz;
%         Cyz_sm = Cyz(1:n,1:m);
%         Cyz_sm_loc = rho_yz(1:n,1:m).*Cyz_sm;
%         kalm = Cyz_sm/(Cz+CV(1:m,1:m));
%         kalm_loc = Cyz_sm_loc/(Cz_loc+CV(1:m,1:m));
%         
%         figure
%         subplot(3,3,1), imagesc(Cz), title('Cz'), colorbar, colormap(bluewhitered)
%         subplot(3,3,2), imagesc(inv(Cz + CV(1:m,1:m))), title('inv Cz'), colorbar, colormap(bluewhitered)
%         subplot(3,3,3), imagesc(kalm), title('K'), colorbar, colormap(bluewhitered)
%         subplot(3,3,4), imagesc(Cz_loc), title('Cz loc'), colorbar, colormap(bluewhitered)
%         subplot(3,3,5), imagesc(inv(Cz_loc + CV(1:m,1:m))), title('inv Cz loc'), colorbar, colormap(bluewhitered)
%         subplot(3,3,6), imagesc(kalm_loc), title('K loc'), colorbar, colormap(bluewhitered)
%         subplot(3,3,7), imagesc(rho_yz(1:n,1:m)), title('\rho_{yz}'), colorbar
%         subplot(3,3,8), imagesc(rho_zz(1:m,1:m)), title('\rho_{zz}'), colorbar
%         
%         Czz_orig = Czz;
%         Czz(50,1) = 0.3;
%         kalm = Cyz/(Czz+CV);
%         kalm_loc = Cyz_loc/(Czz_loc+CV);
%         figure
%         subplot(3,3,1), imagesc(Czz), title('Cz'), colorbar, colormap(bluewhitered), caxis([0,1])
%         subplot(3,3,2), imagesc(inv(Czz + CV)), title('inv Cz'), colorbar, colormap(bluewhitered)
%         subplot(3,3,3), imagesc(kalm), title('K'), colorbar, colormap(bluewhitered)
%         subplot(3,3,4), imagesc(Czz_loc), title('Cz loc'), colorbar, colormap(bluewhitered), caxis([0,1])
%         subplot(3,3,5), imagesc(inv(Czz_loc + CV)), title('inv Cz loc'), colorbar, colormap(bluewhitered)
%         subplot(3,3,6), imagesc(kalm_loc), title('K loc'), colorbar, colormap(bluewhitered)
%         subplot(3,3,7), imagesc(rho_yz), title('\rho_{yz}'), colorbar
%         subplot(3,3,8), imagesc(rho_zz), title('\rho_{zz}'), colorbar
%         
%         
%         figure
%         subplot(2,2,1)
%         imagesc(Czz), colorbar, title('Czz')
%         subplot(2,2,2)
%         imagesc(Czz_loc), colorbar, title('Czz_loc')
%         subplot(2,2,3)
%         imagesc(inv(Czz+CV)), colorbar, title('inv Czz'),% caxis([-60,100])
%         subplot(2,2,4)
%         imagesc(inv(Czz_loc+CV)), colorbar, title('inv rho*Czz')%, caxis([-60,100])
%         
%         figure
%         subplot(2,2,1), imagesc(rho_zz_orig), title('\rho_{zz} orig'), colorbar
%         subplot(2,2,2), imagesc(rho_zz), title('\rho_{zz}'), colorbar
%         subplot(2,2,3), imagesc(zz_mask), title('zz mask'), colorbar
%         subplot(2,2,4), imagesc(inv(Czz_loc_orig+CV)), title('inv Czz loc orig'), colorbar
%         
%         figure
%         subplot(1,2,1)
%         imagesc(kalm_loc), title('kalman gain (local)'), colorbar
%         subplot(1,2,2)
%         imagesc(kalm), title('kalman gain'), colorbar
%         
%         CV = eye(120)*0.01;
%         kalm_loc = Cyz_loc/(Czz_loc + CV);
        
%         k_mask = H'/(H*H');
%         yz_mask(yz_mask>0) = 1;
%         zz_mask(zz_mask>0) = 1;
%         rho_yz_orig = rho_yz;
%         rho_zz_orig = rho_zz;
%         rho_yz = rho_yz_orig;
%         rho_k = rho_yz/rho_zz;
%         rho_yz(logical(yz_mask)) = 1;
%         rho_zz(logical(zz_mask)) = 1;
% %         [ind1,ind2] = find(full(yz_mask)==1);
% %         rho_yz(ind1,ind2) = 1;
%         Cyz_loc = rho_yz.*Cyz;
%         Czz_loc = rho_zz.*Czz;
    end
%     
%     % compute Kalman gain
%     kalm = Cyz/(Czz+CV);
           
    if opt.loc
%         kalm_loc = kalm.*rho_yz; % nah, this ain't it...
%         kalm_loc = Cyz_loc/(Czz+CV);
        kalm_loc = Cyz_loc/(Czz_loc+CV);
%         kalm_loc = Cyz_loc*pinv(Czz_loc+CV);
        if sum(missing_rows)>0 % account for missing rows
            kalm_loc(:,missing_rows) = 0;
        end    
    else
        try
            kalm = Cyz*pinv(Czz+CV);
        catch
            warning('predicted_meas are nan')
            kalm = zeros(size(Cyz)); 
        end
        if sum(missing_rows)>0 % account for missing rows
            kalm(:,missing_rows) = 0;
        end            
    end
    
%     figure
%     imagesc(kalm_loc), title('K case'),colorbar
%     figure
%     imagesc(Cyz_loc), title('Cyz case'),colorbar
%     figure
%     imagesc(Czz_loc), title('Czz case'),colorbar
    
    % ^this method is GOOD because you don't need to store the covariance matrix P!!!
    
    plotflag = 0;
    if plotflag
        check_EnKF_matrices(Cyz, Czz, kalm)
        check_EnKF_matrices(Cyz_loc, Czz_loc, kalm_loc)
    end
%     
%     figure
%     subplot(2,2,1), imagesc(Cyz), title('C_{yz}'),colorbar
%     subplot(2,2,2), imagesc(Czz), title('C_{zz}'),colorbar
%     subplot(2,2,3), imagesc(kalm), title('K'),colorbar
%     subplot(2,2,4), imagesc(inv(Czz+CV)), title('inv(C_{zz}+C_v)'),colorbar, colormap(bluewhitered)
    
    % Measurement error realizations needed for update
    v=mvnrnd(zeros(1,size(Z_ENS,1)),CV,Nreps)';
    
    % Additive update
    
    if sum(missing_rows)==0
        ZMEAS(isnan(ZMEAS)) = 0; % so the NaNs don't get in the way...
        
        Y_CHANGE = Y_UPDATE;
        for k=1:Nreps
            if opt.loc
                Y_UPDATE(:,k)=Y_ENS(:,k)+kalm_loc*(ZMEAS+v(:,k)-Z_ENS(:,k));
                Y_CHANGE(:,k) = kalm_loc*(ZMEAS+v(:,k)-Z_ENS(:,k));
            else
                Y_UPDATE(:,k)=Y_ENS(:,k)+kalm*(ZMEAS+v(:,k)-Z_ENS(:,k));
            end
        end
    
    else  
        
        Y_CHANGE = Y_UPDATE;
        for k=1:Nreps
            if opt.loc
                Y_UPDATE(:,k)=Y_ENS(:,k)+kalm_loc*(ZMEAS+v(:,k)-Z_ENS(:,k));
                Y_CHANGE(:,k) = kalm_loc*(ZMEAS+v(:,k)-Z_ENS(:,k));
            else
                ZMEAS(missing_rows,:) = 0; % handling missing data
                Y_UPDATE(:,k)=Y_ENS(:,k)+kalm*(ZMEAS+v(:,k)-Z_ENS(:,k));
                Y_UPDATE(missing_rows,k)=Y_ENS(missing_rows,k);
            end
        end
        
    end
    
    % compare prior and posterior to measurements
    if plotflag
    m = 95;
    n = 3681;
    s = 17;
    gg = 1;
    prior_runoff_ensemble = unfliparoom(exp(Y_ENS), n, s);
    prior_mean_discharge = state_model_dumb(mean(prior_runoff_ensemble,3), opt.HH);
    post_runoff_ensemble = unfliparoom(exp(Y_UPDATE), n, s);
    post_mean_discharge = state_model_dumb(mean(post_runoff_ensemble,3), opt.HH);
    discharge_meas = unfliparoo(exp(ZMEAS), m, s);
    
    for gg=1:10
        figure
        plot(prior_mean_discharge(:,gg), 'blue', 'linewidth', 2);
        hold on
        plot(post_mean_discharge(:,gg), 'red', 'linewidth', 2);
        plot(discharge_meas(:,gg), 'black', 'linewidth', 2);
        legend('prior','posterior','true')
        title(num2str(gg))
    end
    end
    
return