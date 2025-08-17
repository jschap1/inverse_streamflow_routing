% Plot ACF
% Based on acf.m function from the MATLAB FEX
%
% INPUTS
% ta = sample ACF
% N = sample size used to calculate ta
% p = maximum lag

function ax = plotacf(ta, p, N)

% bar(1:p, [ta])
bar(0:p, [1;ta])

% Plot rejection region lines for test of individual autocorrelations
% H_0: rho(tau) = 0 at alpha=.05
line([0 p+.5], (1.96)*(1/sqrt(N))*ones(1,2))
line([0 p+.5], (-1.96)*(1/sqrt(N))*ones(1,2))

% Some figure properties
line_hi = (1.96)*(1/sqrt(N))+.05;
line_lo = -(1.96)*(1/sqrt(N))-.05;
bar_hi = 1+.05 ;
% bar_hi = max(ta)+.05 ;
bar_lo = -max(ta)-.05 ;
if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
    axis([0 p+.60 line_lo line_hi])
else
    axis([0 p+.60 bar_lo bar_hi])
end
title({' ','Sample Autocorrelations',' '})
xlabel('Lag Length')
set(gca,'YTick',[-1:.20:1])
% set number of lag labels shown
if (N<28 & N>4)
    set(gca,'XTick',floor(linspace(1,p,4)))
elseif (p>=28)
    set(gca,'XTick',floor(linspace(1,p,8)))
end
set(gca,'TickLength',[0 0])

end
