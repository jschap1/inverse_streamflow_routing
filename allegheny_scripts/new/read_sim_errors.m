% Read in simulated runoff error from CoSMoS
%
% 8/1/2023 JRS

% Allegheny basin is 20x21

cd('/Volumes/HD3/ISR/inverse_streamflow_routing')

M = 2000;
nt = 120;
dim = 21;
wkdir = './allegheny_data/errors/m1a1L40T5_LM';
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

