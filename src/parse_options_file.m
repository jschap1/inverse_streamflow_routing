% Parse options file
%
% Reads a text file with run-time options into MATLAB
% Options must be numeric

function opt = parse_options_file(optionsfile)

%% Read options

fmt = '%s%s';
opt = {'Delimiter',' ', 'CommentStyle', '#'};
[fid,msg]=fopen(optionsfile,'rt');
assert(fid>=3,msg);
options = textscan(fid,fmt,opt{:});
fclose(fid);

opt = struct();
for i=1:length(options{1})
    disp(['Setting ' options{1}{i} '=' options{2}{i}])
    eval(['opt.' options{1}{i} '=' options{2}{i} ';'])
end

return