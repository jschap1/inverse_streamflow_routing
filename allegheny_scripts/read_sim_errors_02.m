% Read in simulated runoff error from CoSMoS
%
% 8/1/2023 JRS

% Allegheny basin is 20x21

M = 200;
nt = 365;
dim = 21;
wkdir = '/Volumes/HD3/ISR/inverse_streamflow_routing/allegheny_data/errs/m1s2L40T5_lognorm';
fnames = dir(fullfile(wkdir, '*sim*.txt'));

tic
alldata = zeros(20, dim, nt, M);
for m=1:M
    A = dlmread(fullfile(wkdir, fnames(m).name));
%     A = readmatrix(fullfile(wkdir, fnames(m).name));
    B = reshape(A(1:nt,:), nt, dim, dim);
    C = permute(B, [2,3,1]);
    D = C(1:20,:,:); % basin area is rect, not square
    alldata(:,:,:,m) = D;
    toc
end

save(fullfile(wkdir, ['simulated_runoff_errors_' num2str(M) '.mat']), 'alldata', '-v7.3')

toc

