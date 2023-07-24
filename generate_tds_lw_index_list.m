function generate_tds_lw_index_list(type, version, fdir)
% generate_lw_index_list(Type: sbm1, sbm2, tswf or rswf, version: 0-99, fdir: folder to save the list) 
% This function generates a readable .txt list with the dates and indexes
% from the database created by get_tds_lw_indexes.



% Parsing cdf type
if ~exist('type','var') || isempty(type)
    column = 3; % default tswf
else
    types = ["sbm1", "sbm2", "tswf", "rswf"];
    if ~any(type == types)
        disp('Wrong type, use "sbm1", "sbm2", "tswf" or "rswf".')
        return
    end
end
switch type
    case "sbm2"
        column = 2;
    case "tswf"
        column = 3;
    case "rswf"
        column = 4;
    case "sbm1"
        column = 5;
end

% Parsing version
if ~exist('version', 'var')
    fprintf('Please specify the version of the list, saving as version -1\n')
    version = -1;
end

% Parsing file directory
if ~exist('fdir', 'var')
    fprintf('Output directory not specified, setting fdir = %s\n', pwd)
    fdir = pwd;
end
if ~exist(fdir, 'dir')
    fprintf('Output directory does not exist, creating folder.\n')
    mkdir(fdir)
end
fdir = [fdir filesep 'tds_lw_indexes_' type '_V' sprintf('%02i', version) '.txt'];



% Loading from database
load('tds_lw_indexes.mat','lw_indexes'); 
fileID = fopen(fdir,'w');
nrows = size(lw_indexes,1);
for row = 1:nrows
    if ~isempty(lw_indexes{row, column})
        indexes = sprintf('%i,',lw_indexes{row,column});
        fprintf(fileID,'%10s %s\n',datestr(lw_indexes{row,1},'yyyy-mm-dd'),indexes(1:end-1));
    end
end
fclose(fileID);



end