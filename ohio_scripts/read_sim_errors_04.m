% Read in simulated runoff error from CoSMoS
%
% 4/28/2023 JRS

% basin is 72x104

M = 3;
nt = 365;
dim = 104;
wkdir = './ohio_data/sample_cosmos_errors';
fnames = dir(fullfile(wkdir, '*sim*.txt'));

tic
alldata = zeros(72, dim, nt, M);
for m=1:M
    A = dlmread(fullfile(wkdir, fnames(m).name));
%     A = readmatrix(fullfile(wkdir, fnames(m).name));
    B = reshape(A(1:nt,:), nt, dim, dim);
    C = permute(B, [2,3,1]);
    D = C(1:72,:,:); % basin area is rect, not square
    alldata(:,:,:,m) = D;
    toc
end

save(fullfile(wkdir, ['simulated_runoff_errors_' num2str(M) '.mat']), 'alldata', '-v7.3')

toc

