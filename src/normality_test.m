% Using kolmogorov-smirnov test for normality

function [isnorm, p] = normality_test(x)

    siglevel = 10; % percent

    dat = x(:); 
    s = (dat - mean(dat))/std(dat);
    
%     figure,histogram(s)
    
    [h,p,ks, cv] = kstest(s, 'Alpha', siglevel/100);
    
    if h==0
        % cannot reject null hypothesis of normdist
        isnorm = 1;
    elseif h==1
        % reject null H0 of normalist
        isnorm = 0;
    end
    
return