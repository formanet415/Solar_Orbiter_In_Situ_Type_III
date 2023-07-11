function [ttpk, f2] = t3_auto_fit_electron_vel(ep0, tlen, species, opts, justfit)

%T3_AUTO_FIT_ELECTRON_VEL does automatic fitting of the electron beam speed
%   from EPD data.
% code from solo_panel_epd_step_main/rates_spectrum

% also saves a plot of the fit
%
% Options (opts)
%  opts.max_freq      - maximum frequency in kHz [def = 80]
%  opts.show_xlabel   - if true, show label on X axis [default = 1]
%  opts.show_colorbar - if true, show colorbar [default = 0]
%  opts.time_seconds  - if true, time axis is labeled in seconds
%                       otherwise datetick is used. [def = 0]

ttpk = [];
f2= [];
% default parameters

show_xlabel = 1;
show_title = 1;
time_datenum = 1;
opts_sp.colorax = 0;
opts_sp.log_y = 1;
threshold = 1;

if ~exist('justfit','var') || isempty(justfit)
    justfit = 0;
end
if ~exist('species','var') || isempty(species)
    species = 'electrons';
end

if exist('opts','var') && ~isempty(opts)
    if isfield(opts,'colorax')
        opts_sp.colorax = opts.colorax;
    end
    if isfield(opts,'justfit')
        justfit = opts.justfit;
    end
    if isfield(opts,'show_xlabel')
        show_xlabel = opts.show_xlabel;
    end
    if isfield(opts,'max_energy')
        opts_sp.max_freq = opts.max_energy;
    end
    if isfield(opts,'flux_threshold')
        threshold = opts.flux_threshold;
    end
    if isfield(opts,'show_xticks')
        opts_sp.show_xticks = opts.show_xticks;
    end
    if isfield(opts,'show_title')
        show_title = opts.show_title;
    end
    if isfield(opts,'show_colorbar')
        if (0 == opts.show_colorbar)
            opts_sp.colorax = 0;
        elseif (0 == opts_sp.colorax)
            opts_sp.colorax = 1;
        end
    end
    if isfield(opts,'time_seconds')
        time_datenum = ~opts.time_seconds;
    end
end


if ep0<datenum(2021,10,22)
    [ep, int_flux, mag_flux, energ] = caadb_get_solo_epd_step_rates(ep0, tlen);
elseif ep0>datenum(2021,10,22)
    [ep, ~,int_flux, ~, mag_flux, energ] = caadb_get_solo_epd_step_main(ep0, tlen);
end

if (isempty(ep))
    plot(0);
    title('No EPD/STEP data for this event');
    return;
end

switch (species)
    case 'electrons'
        spec = int_flux - mag_flux;
        titlabel = 'Electron flux (Int-Mag)';
        axlabel = 'energy [keV]';
    case 'ions'
        spec = mag_flux;
        titlabel = 'Ion flux (Mag)';
        axlabel = 'energy [keV/n]';
    case 'total'
        spec = int_flux;
        titlabel = 'Total flux (Int)';
        axlabel = 'energy [keV]';
    otherwise
        fprintf(1,'Unknown species %s\n', species);
        return;
end

spec(spec <= threshold) = nan;


spn = normalize(spec,2);
fpspn = movmean(spn,5,2);
dims = size(fpspn);
ttpk = [];
dtpk = [];
for i = 1:dims(2)
    [pks, loc] = findpeaks(fpspn(:,i),'MinPeakDistance',7,'MinPeakProminence',1.7);
    peaks(i,loc) = pks;
    ttpk = [ttpk ones(length(pks),1)'*ep(i)];
    dtpk = [dtpk energ(loc)'*1e3];
end
if justfit==0
    figure(1)
    plot_spectrogram(spn, energ*1e3, ep, opts_sp);
    %figure(2)
    %imagesc(peaks')
    %set(gca,'YDir','normal')
    %figure(3)
    %scatter(ttpk,dtpk)
    %datetick()
end
try
    f1 = fit(ttpk',log(dtpk)','poly1');
catch ME
    disp('fitting failed')
    return
end
%hold on
%plot(ttpk,exp(f1(ttpk)))
%hold off
x=dtpk'-exp(f1(ttpk));
filter = abs(x)-std(x)<0;
if sum(filter)<10
    return
end
f2 = fit(ttpk(filter)',log(dtpk(filter))','poly1','Weights',exp(linspace(5,1,sum(filter))));
%figure(1)
if f2.p1>0
    disp('fitting failed, p1 positive')
    return
end
hold on
scatter(ttpk,dtpk,'black')
plot(ttpk,exp(f2(ttpk)),'b','LineWidth',5);
if justfit == 0
    xlim([ep0 ep0 + tlen/86400]);

    set(gca,'FontSize',12);
    ylabel(axlabel,'FontSize',12);
    if (show_title)
        title(sprintf('EPD STEP-MAIN omnidirectional %s %s', titlabel, datestr(ep0,'YYYY-mm-dd HH:MM:SS.FFF')),'FontSize',12);
    end
    if ~show_xlabel
        xlabel([]);
    end
    datetick()
    hold off
end
end

