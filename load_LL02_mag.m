function llmg = load_LL02_mag(date)
%LOAD_LL02_MAG Loads low latency MAG data from the .cdf files

% finding file
if isa(date, 'double')
    date = datetime(date,'ConvertFrom','datenum');
end
[y, m, d] = ymd(date);

homedir = getenv('OKF');

rocstr = sprintf('solo_LL02_mag_%04i%02i%02iT', y , m, d);
fdir = [homedir 'mag' filesep 'LL02' filesep sprintf('%04i',y) filesep]; 
files = dir(fdir);
files = struct2cell(files); files = files(1,:);

for i = 1:length(files)
    file = files{i};
    if contains(file, rocstr)
        fname = file;
        fprintf('found LL02 mag file: "%s"\n', fname)
    end
end

if ~exist('fname') || isempty(fname)
    fprintf('LL02 data on %04i%02i%02i not found\n', y , m, d)
    llmg = [];
    return
end

% loading file

[cstr, cinfo] = spdfcdfread([fdir filesep fname]);

llmg.epoch = cstr{1};
llmg.data_rtn = cstr{2};
llmg.rtn_label = cstr{2};
llmg.data_srf = cstr{6};
llmg.srf_label = cstr{7};
end

