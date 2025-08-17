function [sub1runoff, sub2runoff] = calc_med_runoff_nse_for_localiz_fig(sub1runoff, sub2runoff, gi)

n1 = 233-50;
n2 = 50;

sub1runoff.nse_prior = zeros(n1,1);
for j=1:n1
    sub1runoff.nse_prior(j) = myNSE(sub1runoff.true(j,gi), sub1runoff.prior(j,gi));
    sub1runoff.nse_us(j) = myNSE(sub1runoff.true(j,gi)', sub1runoff.us(j,gi)');
    sub1runoff.nse_ds(j) = myNSE(sub1runoff.true(j,gi), sub1runoff.ds(j,gi));
    sub1runoff.nse_both(j) = myNSE(sub1runoff.true(j,gi), sub1runoff.both(j,gi));
end

% downstream sub-basin median, cell-wise runoff improvement
median(sub1runoff.nse_us') - median(sub1runoff.nse_prior') % 0.11
median(sub1runoff.nse_ds') - median(sub1runoff.nse_prior') % 0.28
median(sub1runoff.nse_both') - median(sub1runoff.nse_prior') % 0.28

sub2runoff.nse_prior = zeros(n2,1);
for j=1:n2
    sub2runoff.nse_prior(j) = myNSE(sub2runoff.true(j,gi), sub2runoff.prior(j,gi));
    sub2runoff.nse_us(j) = myNSE(sub2runoff.true(j,gi), sub2runoff.us(j,gi));
    sub2runoff.nse_ds(j) = myNSE(sub2runoff.true(j,gi), sub2runoff.ds(j,gi));
    sub2runoff.nse_both(j) = myNSE(sub2runoff.true(j,gi), sub2runoff.both(j,gi));
end

% upstream sub-basin median, cell-wise runoff improvement
median(sub2runoff.nse_us') - median(sub2runoff.nse_prior') % 0.14
median(sub2runoff.nse_ds') - median(sub2runoff.nse_prior') % 0.15
median(sub2runoff.nse_both') - median(sub2runoff.nse_prior') % 0.14

return