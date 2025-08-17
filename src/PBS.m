% Uses the particle batch smoother to perform the ISR update
%
% 9/22/2023 JRS

function ypost_mean = PBS(Y_ENS,Z_ENS,ZMEAS,Cv, missing_rows, opt)

   if all(missing_rows==1)
        ypost_mean = Y_ENS;
        return
   end
   
   % Remove missing measurements
    Z_ENS(missing_rows,:) = [];
    ZMEAS(missing_rows,:) = [];
   
    % Get size of ensemble
    nparticles = size(Y_ENS,2);
      
    % Calculate likelihood
    %pred_meas = squeeze(zz(1,end,:));
    %actual_meas = z(jmeas);
    
    % likelihood function
    % number of observations in batch
    Nobs=length(ZMEAS);
    % Covariance matrix for batch
    C_VV=diag(Cv*ones(Nobs,1));
    for k=1:nparticles
        L(k) = 1/(2*pi)^Nobs/det(C_VV)*...
            exp(-1/2*(ZMEAS(:)-Z_ENS(:,k))'*inv(C_VV)*...
                                        (ZMEAS(:)-Z_ENS(:,k)));
    end

    % Normalize likelihood of each particle
    Lsum = sum(L);
    for ll = 1:1:nparticles
        L(ll) = L(ll) / Lsum;
    end
    posterior_weights=L;
    
    % Compute posterior statistics
    ypost_mean = Y_ENS*posterior_weights'; 
    tmp = Y_ENS.^2;
    ypost_std = (posterior_weights*tmp')' - ypost_mean.^2;
     

return