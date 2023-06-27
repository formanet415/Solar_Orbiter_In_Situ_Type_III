function tds_analyze_type3(filename, type, opts)
% t3_analyze_sbm2 This function performs analysis on the wave polarisation
% during a type III event recorded in SBM2 mode by SolO
%   Use:
%    Important opts settings:
%       mode - used for Langmuir wave detection, see GET_TDS_LW_INDEXES.
%
%    Other opts settings:
%       hrs - how many hours to include in the plot/analysis
%
%   Requirements: MAG data


% List of figures:
% - TNR + EPD + wave polarisation + angle between B and the Y-Z plane
% - wave polarisation vs. beam energy (color ~ quality ~ angle B vs. Y-Z)


% Loading data
caa_data_paths
cdf = tdscdf_load_l2_surv_tswf(filename);
ep0 = cdf.epoch(1);
%rswf = tdscdf_load_l2_surv_rswf(ep0, 1);


% Setting parameters
if exist("opts") %#ok<EXIST>
    if isfield(opts, 'hrs')
        tlen = opts.hrs*3600;
    else
        tlen = 2.5*3600;
    end

    if isfield(opts, 'tnr_comp')
        tnr_comp = opts.tnr_comp;
    else
        tnr_comp = 2;
    end

    if isfield(opts, 'mode')
        mode = opts.mode;
    else
        mode = 0;
    end

    if isfield(opts, 'tnr_time')
        tnr_time = opts.tnr_time;
    else
        tnr_time = 0;
    end

    if isfield(opts, 'epd_time')
        epd_time = opts.epd_time;
    else
        epd_time = 0;
    end
else 
    mode = 0;
    tnr_time = 0;
    epd_time = 0;
    tnr_comp = 2;
end

if ~exist('type','var') || isempty(type)
    type = "sbm2";
else
    types = ["sbm2", "tswf", "rswf"];
    if ~any(type == types)
        disp('Wrong type, use "sbm2", "tswf" or "rswf".')
    end
end


% Detect Langmuir waves
indexes = get_tds_lw_indexes(cdf,type,mode);
if isempty(indexes)
    disp('No Langmuir waves detected in this file.')
end

% Determine the Langmuir wave polarization
f = zeros(size(indexes));
quality = zeros(size(indexes));
fep = cdf.epoch(indexes);
t_avg = 0.5/86400;
for i = 1:length(indexes)
    ind = indexes(i);
    fs = cdf.samp_rate(ind);
    % Convert waveform to SRF
    wf(1:2,:) = convert_to_SRF(cdf,ind);
    
    % Loading MAG data
    [~, b, ~, ~, ~]=caadb_get_solo_mag(fep(i)-t_avg/2,t_avg*86400,'srf');
    if isempty(b)
        fprintf('Mag data missing, looking for low latency data\n')
        [~, b, ~, ~] = caadb_get_solo_mag_LL02(fep(i)-t_avg/2,t_avg*86400,'srf');
        if isempty(b)
            fprintf('NO MAG DATA FOUND\n')
            continue
        end
    end

    % Project B onto the Y-Z plane
    b = mean(b,2);
    b_yz = b(2:3);
    qual = vecnorm(b_yz)/vecnorm(b);
    quality(i) = qual;

    % FAC transformation
    e_par = b_yz/vecnorm(b_yz);
    e_ort = [e_par(2); -e_par(1)];

    w_par = e_par'*wf/qual;
    w_ort = e_ort'*wf;
    
    % find Langmuir wave peaks for bandpower
    s_par = fft(w_par); s_ort = fft(w_ort);
    n = length(w_par);                                  % number of samples
    freq = (0:n-1)*(fs/n);                              % frequency range
    mask = ((freq>10e3) + (freq<70e3) == 2);
    s_par = abs(s_par(mask)).^2/n; s_ort = abs(s_ort(mask)).^2/n;   % power of the DFT
    freq = freq(mask);
    [~,locs_per,~,~] = findpeaks(s_par,'SortStr','descend');  
    [~,locs_ort,~,~] = findpeaks(s_ort,'SortStr','descend'); 
    lwfq = mean([freq(locs_per(1)),freq(locs_ort(1))]);

    % calculate energy of the Langmuir wave components
    p_par = bandpower(s_par,fs,[lwfq-2.5e3 lwfq+2.5e3]);
    p_ort = bandpower(s_ort,fs,[lwfq-2.5e3 lwfq+2.5e3]);
    f(i) = p_ort/(p_par+p_ort);
end

if tnr_time == 0
    % Selecting TNR_time
    ep1 = cdf.epoch(end);
    tt = 86400*(ep1-ep0);
    solo_panel_tnr_spectrum(ep0-1/24,tt+1*3600,tnr_comp);
    disp('Select the beginning of the radio emission')
    [tnr_time,~] = ginput(1);
    fprintf('TNR time set to %s\n', datestr(tnr_time))
end

if epd_time == 0
    % Selecting EPD_time
    if ep0<datenum(2021,10,22)
        solo_panel_epd_step_rates_spectrum(tnr_time-0.25/24,3*3600,'electrons')
    elseif ep0>datenum(2021,10,22)
        solo_panel_epd_step_main_spectrum(tnr_time-0.25/24,3*3600,'electrons')
    end
    disp('Select the beginning of the electron beam')
    [epd_time,~] = ginput(1);
    fprintf('EPD time set to %s\n', datestr(epd_time))
end

justfit.justfit = 1;
[ttpks, fit] = t3_auto_fit_electron_vel(epd_time-1/48,3600*2,'electrons',justfit);
fitep0 = epd_time;
fitep1 = ttpks(end);
energs = zeros(size(f));
for i = 1:length(fep)
    if ((fep(i)>fitep0)+(fep(i)<fitep1)) == 2
        energs(i) = exp(fit(fep(i)));
    else
        energs(i) = nan;
    end
end



% Subplot setup
tsubplot = @(m,n,p) subtightplot (m,n,p,[0.06 0.06], [0.05 0.07], [0.07 0.07]);
figure(1); clf;
k = 1; nfig = 1;






end