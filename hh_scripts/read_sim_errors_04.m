% Read in simulated runoff error from CoSMos
%
% 1/10/2023 JRS

M = 500;
nt = 722;
dim = 7;
% wkdir = '/hdd/ISR/02_Basins/HH2/Data/prior_err_L1_S24_mu2_sigma_2/';

% fname <- paste0("L", L_vals[i], "_T", T_vals[j], "_mu", mu_vals[k], "_sigma", sigma_vals[k])
% dirname <-paste0("/hdd/ISR/02_Basins/HH2/Data/runoff_priors/", fname)
        

L_vals = [0.01, 1.57, 7, 14, 21, 28, 56, 83, 538];
T_vals = [0.01, 1.5, 1.875, 3, 4.5, 6, 12, 18, 117];
COV = 1;
mu_vals = [0.5, 1, 2];
sigma_vals = mu_vals/COV;

for i=1:length(L_vals)
    for j=1:length(T_vals)
        for k=1:length(mu_vals)
            
dirname = ['L', num2str(L_vals(i)), '_T', num2str(T_vals(j)), ...
    '_mu', num2str(mu_vals(k)), '_sigma', num2str(sigma_vals(k))];
wkdir_orig = '/hdd/ISR/02_Basins/HH2/Data/runoff_priors/';
wkdir = fullfile(wkdir_orig, dirname);

% wkdir_orig = '/hdd/ISR/02_Basins/HH2/sensitivity_to_prior/varyingSIGMA/';
% dirnames = dir([wkdir_orig 'sigma*']);
    
% wkdir = fullfile(wkdir_orig, dirnames(testnum).name);

fnames = dir(fullfile(wkdir, '*sim*.txt'));

if exist(fullfile(wkdir, 'simulated_runoff_errors_500.mat'))>0
    continue
end

tic
alldata = zeros(dim-1, dim, nt, M);
for m=1:M
    A = readmatrix(fullfile(wkdir, fnames(m).name));
    B = reshape(A(1:nt,:), nt, dim, dim);
    C = permute(B, [2,3,1]);
    D = C(1:6,:,:); % basin area is rect, not square
    alldata(:,:,:,m) = D;
    
    % Delete unneeded files to save storage space
    delete(fullfile(wkdir, fnames(m).name))
    
end

save(fullfile(wkdir, 'simulated_runoff_errors_500.mat'), 'alldata')
disp('Finished reading errors in directory:')
disp(dirname)

toc

        end
    end
end

% Delete unneeded files
