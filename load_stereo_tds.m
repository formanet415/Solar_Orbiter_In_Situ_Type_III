function ster = load_stereo_tds(day,sc)
% LOAD_STEREO_TDS Loads stereo data from the binary files on jaruska
%   ster = load_stereo_tds(day,sc) loads data to structure.
%       day is of datenum or datetime format
%       sc is either "A", "a", "B" or "b" 
%       


 
% Modify path below for own use
% jaruskadir = 'E:\Dokumenty\OneDrive - Univerzita Karlova\OKF\jaruska_offline';
load('jaruskadir.mat') % added to gitignore (contains a single string with path to jaruska)


% when spacecraft is not specified TBD (not practical)
%if ~exist('sc') || isempty(sc)
%    mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
%    a = load_stereo_tds(day,'a');
%    b = load_stereo_tds(day,'b');
%    ster = mergestructs(a,b);      struct fields would need to be
%                                   different for sta and stb
%end

% finding file
if isa(day, 'double')
    day = datetime(day,'ConvertFrom','datenum');
end
[y, m, d] = ymd(day);
sc = lower(sc);
fname = [jaruskadir filesep 'mnt' filesep 'raid' filesep 'stereo' filesep 'swaves_tds' filesep num2str(y) filesep sprintf('stereo_tds_%s_%04i-%02i-%02i.bin', sc, y , m, d)];


fileID = fopen(fname, 'rb');
if fileID == -1
    fprintf('file not found: %s', fname);
    ster = [];
    return
end


nsnap = fread(fileID, 1, 'int32');
sn_len = fread(fileID, 1, 'int32');

yr(1:nsnap) = fread(fileID, nsnap, 'int32');
doy(1:nsnap) = fread(fileID, nsnap, 'int32');
hr(1:nsnap) = fread(fileID, nsnap, 'double');
antenna(1:nsnap) = fread(fileID, nsnap, 'int32');
dt(1:nsnap) = fread(fileID, nsnap, 'double');
wlen(1:nsnap) = fread(fileID, nsnap, 'int32');
filter(1:nsnap) = fread(fileID, nsnap, 'int32');
channel(1:nsnap) = fread(fileID, nsnap, 'int32');
data(1:sn_len,1:nsnap) = fread(fileID, [sn_len, nsnap], 'float');

ster.ep(1:nsnap/4) = datenum(yr(1:4:nsnap),1,1) + doy(1:4:nsnap) - 1 + hr(1:4:nsnap)/24;
ster.n_samps(1:nsnap/4) = wlen(1:4:nsnap);
for i = 1:4    
    ster.data(1:sn_len,i,1:nsnap/4) = data(:,i:4:nsnap);
    ster.channel(i,1:nsnap/4) = channel(i:4:nsnap);
    ster.antenna(i,1:nsnap/4) = antenna(i:4:nsnap);
end
ster.samp_rate(1:nsnap/4) = 1e3./dt(1:4:nsnap);
ster.mod(1:nsnap/4) = 60.*hr(1:4:nsnap);
ster.sc = sc;

end

