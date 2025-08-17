% Plot runoff map snapshots
%
% INPUTS
% t = timesteps where we are plotting snapshots
% cbnds = colorbar bounds for each snapshot (nx2 matrix)
%
% Also computes median NSE

function [NSE] = plot_runoff_map_snapshots(t, cbnds, basin, PW13_runoff, Y20_runoff, ENS_runoff, tmpa, nldas)

figure
nt = length(t);
subi =1;

n = size(Y20_runoff,1);

for i=1:nt

%     prior_runoff_map = make_map(basin, tmpa_runoff_prior(t(i), :));
%      post_runoff_map_PW13 = make_map(basin, PW13_runoff(t(i),:));
     post_runoff_map_Y20 = make_map(basin, Y20_runoff(t(i),:));
     post_runoff_map_ENS = make_map(basin, ENS_runoff(t(i),:));
    
%      post_runoff_map_PW13(post_runoff_map_PW13<0) = -999;
     post_runoff_map_Y20(post_runoff_map_Y20<0) = -999;
     post_runoff_map_ENS(post_runoff_map_ENS<0) = -999;
     
    % Prior
    subplot(nt,4,subi)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['TMPA Prior (day ' num2str(t(i)) ')'])
    caxis([cbnds(i,:)])
    
    % Posterior (PW13)
%     subplot(5,nt,nt+i)
%     plotraster(basin.lonv, basin.latv, post_runoff_map_PW13, ['PW13 Posterior (day ' num2str(t(i)) ')'])
%     caxis([cbnds(i,:)])
    
    % Posterior (Y20)
    subi = subi+1;
    subplot(nt,4,subi)
    plotraster(basin.lonv, basin.latv, post_runoff_map_Y20, ['Y20 Posterior (day ' num2str(t(i)) ')'])
    caxis([cbnds(i,:)])
    
    % Posterior (Ensemble)
    subi = subi+1;
    subplot(nt,4,subi)
    plotraster(basin.lonv, basin.latv, post_runoff_map_ENS, ['ENS Posterior (day ' num2str(t(i)) ')'])
    caxis([cbnds(i,:)])
    
    % Truth
    subi = subi+1;
    subplot(nt,4,subi)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['NLDAS Truth (day ' num2str(t(i)) ')'])
    caxis([cbnds(i,:)])
    
    subi=subi+1;
    
    NSE.tmpa(i) = median(myNSE(basin.true_runoff(t(i),:)', basin.prior_runoff(t(i),:)'));
    NSE.y20(i) = median(myNSE(basin.true_runoff(t(i),:)', Y20_runoff(t(i),:)'));
    NSE.ens(i) = median(myNSE(basin.true_runoff(t(i),:)', ENS_runoff(t(i),:)'));
    
end
colormap cool

return

