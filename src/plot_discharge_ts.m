function kge = plot_discharge_ts(post_discharge, true_discharge)

% Should also compute NSE and plot it on the maps
% Also should give the ability to choose which gauges to plot
% Also should allow user to pass in a time vector

[nt,m]=size(post_discharge);
lw=2;
fs=18;

% figure
if m==1
    plot(true_discharge, 'k','linewidth', lw)
    hold on
    plot(post_discharge, 'r', 'linewidth', lw)
    title(['Discharge'])
    legend('Posterior', 'Truth')
    xlabel('DOY')
    ylabel('Discharge (mm/day)')
    set(gca, 'fontsize', fs)  
elseif m<=3
    for i=1:m
        subplot(1,m,i)
        plot(post_discharge(:,i), 'r','linewidth', lw)
        hold on
        plot(true_discharge(:,i), 'k','linewidth', lw)
        title(['Discharge at gauge ' num2str(i)])
        legend('Posterior', 'Truth')
        xlabel('DOY')
        ylabel('Discharge (mm/day)')
        set(gca, 'fontsize', fs)        
    end
elseif m>3
    msub = randperm(m,3);
    for i=1:3
        subplot(1,3,i)
        plot(post_discharge(:,i), 'r', 'linewidth', lw)
        hold on
        plot(true_discharge(:,i), 'k', 'linewidth', lw)
        title(['Discharge at gauge ' num2str(msub(i))])
        legend('Posterior', 'Truth')
        xlabel('DOY')
        ylabel('Discharge (mm/day)')
        set(gca, 'fontsize', fs)        
    end
end

kge = zeros(m,1);
for i=1:m
    kge(i) = myKGE(true_discharge(:,i), post_discharge(:,i));
end

return