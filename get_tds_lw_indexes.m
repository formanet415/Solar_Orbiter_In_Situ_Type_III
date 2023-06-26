function [indexes] = get_tds_lw_indexes(cdf, type, mode)
%GET_TDS_LW_INDEXES looks for langmuir waves in the 
%   This function makes a spectrum of each waveform and looks for peaks
%   characteristic for langmuir waves. 
%
%
%  Input
%   cdf - the tds data in a structure
%
%   type - 'sbm2', 'tswf', 'rpw' (used to categorize data in the database),
%    if left empty, mode 2 is triggered
%
%   mode - There are three modes.
%    mode = 0 - default mode (looks for the indexes in the database (TBD)
%    and if they are not there it performs an automatic detection
%    mode = 1 - manual mode where each pick is checked and the indexes can
%    then be written into the database
%    mode = 2 - overrides the database and performs an automatic detection

if ~exist('type','var') || isempty(type)
    column = 0;
else
    types = ["sbm2", "tswf", "rswf"];
    if ~any(type == types)
        disp('Wrong type, use "sbm2", "tswf" or "rswf".')
    end
end
switch type
    case "sbm2"
        column = 2;
    case "tswf"
        column = 3;
    case "rswf"
        column = 4;
end



if ~exist('mode','var') || isempty(mode)
    mode = 0;
end
indexes = [];

% Load from database
load('tds_lw_indexes.mat'); %#ok<LOAD>
ep0 = cdf.epoch(1);
if mode == 0
    [~,row] = find([lw_indexes{:}] == ep0); %#ok<NODEF> 
    if isempty(row)
        disp('Indexes not found in database, performing automatic detection.')
        mode = 2;
    else
        indexes = lw_indexes{row,column};
        return
    end
end

tsubplot = @(m,n,p) subtightplot (m,n,p,[0.06 0.06], [0.05 0.07], [0.07 0.07]);

% Automatic detection
wfs = cdf.data;
maxsamps = size(wfs, 2);
nsnaps = size(wfs, 3);
fs = cdf.samp_rate;

 

% Preallocate spectrum array
spectra = zeros(3, maxsamps, nsnaps);

% Compute spectra for each waveform
for i = 1:nsnaps
    for j = 1:3
        % Extract waveform for current channel and waveform index
        wf = wfs(j, 1:cdf.samples_per_ch, i);
        
        % Compute spectrum using FFT
        spec = fft(wf);
        n = length(wf);          % number of samples
        f = (0:n-1)*(fs(i)/n);     % frequency range
        spec = abs(spec).^2/n;    % power of the DFT
        if mode==1
            prom = 0.000005;
        elseif mode == 2
            prom = 0.00005;
        end
        [pks,locs,w,~] = findpeaks(spec,'MinPeakProminence',prom, 'MinPeakDistance',50, 'MinPeakWidth', 1); %#ok<ASGLU> 
        freqs = f(locs);
        mask = (freqs<70e3) + (freqs>10e3) == 2;
        freqs = freqs(mask);
        %disp(w(mask))
        if length(pks(mask))>10 && descendingSortingScore(pks(mask))>0.8 % eliminates strong dust
            break
        end

        % Detction
        if ~isempty(freqs)
            if (mode == 2) && (isempty(indexes) || (indexes(end) ~= i))
                indexes(end+1) = i; %#ok<AGROW> 
            elseif (mode == 1)  && (isempty(indexes) || (indexes(end) ~= i))
                tsubplot(2,1,1)
                plot(wf)
                ylabel('E (V/m)')
                xlabel('samples')
                
                tsubplot(2,1,2)
                plot(f*1e-3,log(spec))
                xlim([0,90])
                xlabel('frequency (kHz)')
                ylabel('power')
                for ii = freqs
                    vertline(ii*1e-3,'red');
                end
                
                % User input to check detection
                fprintf('(%i/%i)\n',i,nsnaps)
                answer = questdlg(sprintf('Is this a Langmuir wave?\nTime: %s',datestr(cdf.epoch(i))), ...
                    'Langmuir wave detection', ...
                    'Yes','No','No');
                % Handle response
                switch answer
                    case 'Yes'
                    fprintf('Langmuir wave on index %i\n',i)
                    indexes(end+1) = i; %#ok<AGROW> 
                case 'No'
                    % do nothing
                end
                
            end
        end
        
        % Store the spectrum in the spectra array
        spectra(j, 1:length(spec), i) = spec;
    end
end


if mode == 1
    [row,~] = find([lw_indexes{:}] == ep0);  
    if ~isempty(row)
        lw_indexes{row, column} = indexes;  %#ok<NASGU> 
    else
        lw_indexes{end+1,1} = ep0;
        lw_indexes{end, column} = indexes;

        save("tds_lw_indexes.mat","lw_indexes");
    end
end

end