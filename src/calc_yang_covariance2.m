% Calculates Y20 covariance
%
% Updated 10/11/2022 JRS to use memory efficiently
% Option to supply rho_x so it doesn't have to be re-calculated

function Cx = calc_yang_covariance2(basin, opt, n, s, x, varargin)

if nargin == 5

    % Calculate SC and TC
   
    % Use assumed values of L and T (to generate Cx for ISR)
    opt.k1 = 1/opt.L; % per km
    opt.k2 = 1/opt.T; % per day
 
    % if necessary can load distmat from elsewhere
%     aa = load('./basin_setup/setup-1-gage.mat');
%     basin.distmat = aa.basin.distmat;
    
    SC = calc_SC(opt.k1, basin.distmat, opt.rho_thres, n, s, opt.RAM);
    save(opt.SC_savename, 'SC', '-v7.3')
    
%     figure
%     imagesc(SC(1:2*n,1:2*n))
%     colorbar
        
%     load('./covariance_matrices/TC_1.mat');
    TC = calc_TC(opt.k2, opt.rho_thres, n, s, opt.RAM);
    save(opt.TC_savename, 'TC', '-v7.3')
    
%     figure
%     imagesc(TC(1:5*n,1:5*n))
%     colorbar

    rho_x = SC.*TC;
    save(opt.rho_savename, 'rho_x', '-v7.3')
    clearvars SC TC      
    
%     figure
%     imagesc(rho_x(1:5*n,1:5*n))
%     colorbar
    
elseif nargin == 6
    
    % Re-use rho_x
    
    rho_x = load(varargin{1});    
    rho_x = rho_x.rho_x;
    
elseif nargin == 7
    
    % Calculating Cx for runoff prior generation

    opt.k1 = 1/opt.tL; % per km
    opt.k2 = 1/opt.tT; % per day
 
    SC = calc_SC(opt.k1, basin.distmat, opt.rho_thres, n, s, opt.RAM);
    save(opt.tSC_savename, 'SC', '-v7.3')
    
    TC = calc_TC(opt.k2, opt.rho_thres, n, s, opt.RAM);
    save(opt.tTC_savename, 'TC', '-v7.3')
    
    rho_x = SC.*TC;
    save(opt.trho_savename, 'rho_x', '-v7.3')
    clearvars SC TC  
    
end

% Calculate phi

% storage requirements
S = 4*length(x)^2/1e9;

if opt.progress
    disp(['Single-precision phi will require ' num2str(S) ' GB of memory']);
end

if S < opt.RAM

    phi = calc_phi(x, s, opt.alpha);
    Cx = phi.*rho_x;

else

    disp('phi does not fit in memory')
    disp(' we will need to calculate Cx without explicitly storing phi')

    % Find indices of nonzero values in Cx
    [i,j] = find(rho_x~=0);
    lin_ind = logical(rho_x(:)~=0);

    % Calculate nonzero values of Cx
    nonzero_vals = opt.alpha^2*rho_x(lin_ind).*x(i).*x(j);

    % Assemble sparse matrix from indices and values
    Cx = sparse(i,j,nonzero_vals, n*(s+1), n*(s+1));

end
    
return

%% Extra code (not used)

% % Generate covariance/correlation plots for first time window
% if opt.plotcov
%     figure
%     subplot(2,2,1)
%     imagesc(SC), title('SC')
%     subplot(2,2,2)
%     imagesc(TC), title('TC')
%     subplot(2,2,3)
%     imagesc(phi), title('phi')
%     subplot(2,2,4)
%     imagesc(Cx), title('Cx')
% end

% Cx(Cx<opt.rho_thres) = 0; % (redundant) localization (already applied..)