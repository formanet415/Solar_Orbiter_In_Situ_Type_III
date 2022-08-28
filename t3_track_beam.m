function [ep, bpathlen, hpathlen] = t3_track_beam(tnr_time, tlen, opts)

%  modified code from solo_panel_tnr_spectrum(ep0, tlen, sc, opts)
%
% t3_track_beam(tnr_time, tlen, component, opts)
% Traces beam path RPW TNR spectrogram
%
%
% Options (opts)
%  opts.max_freq      - maximum frequency in kHz [def = 80]
%  opts.show_xlabel   - if true, show label on X axis [default = 1]
%  opts.show_colorbar - if true, show colorbar [default = 0]
%  opts.time_seconds  - if true, time axis is labeled in seconds
%                       otherwise datetick is used. [def = 0]


% default parameters

show_xlabel = 1;
show_title = 1;
time_datenum = 1;
opts_sp.colorax = 0;
opts_sp.log_y = 1;

opts_sp.regular_time  = 0;
opts_sp.max_bin_width = 1/(24*60);


component = 1;

if exist('opts','var') && ~isempty(opts)
    if isfield(opts,'colorax')        
		opts_sp.colorax = opts.colorax;    
    end
    if isfield(opts,'min_freq')        
		opts_sp.min_freq = opts.min_freq;    
    end	
	if isfield(opts,'show_xlabel')
		show_xlabel = opts.show_xlabel;
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

[ep, spec, fq, axlabel] = cdb_solo_build_thr_spectrogram(tnr_time, tlen, component);

if (isempty(ep))
    plot(0);
    title('No TNR data for this event');
    return;
end

effl = 3.24;
spn = normalize(log(spec(60:end,:)),2);
spn = movmean(spn,10,2);
sfq = fq(60:end);

dims = size(spn);
ttpk = [];
dtpk = [];
for i = 1:dims(1)
    [pks, loc] = findpeaks(spn(i,:),'MinPeakDistance',10,'MinPeakProminence',1);
    %[pks, loc] = max(spn(:,i));
    peaks(i,loc) = pks;
    ttpk = [ttpk mean(ep(loc))];
    dtpk = [dtpk sfq(i)];
end

[~, ~, vel_rtn] = caadb_get_solo_swa_pas_moments(tnr_time,tlen);
tmp = vecnorm(vel_rtn);
sw_vel = mean(tmp);

figure(3)
plot_spectrogram(spec*(effl*effl)^2, fq/1e3, ep, opts_sp);
%xlim([ep0 ep0 + tlen/86400]);

set(gca,'FontSize',12);
ylabel('frequency [kHz]','FontSize',12);
if (show_title)
    title(sprintf('RPW TNR %s %s', axlabel, datestr(tnr_time,'YYYY-mm-dd HH:MM:SS.FFF')),'FontSize',12);
end
if ~show_xlabel
	xlabel([]);
end
hold on
plot(ttpk,dtpk*1e-3,'ro')

hold off

[bd, hd] = density_model_simple(dtpk);
bpath_length = parker_spiral_length(bd/1.496e+11,sw_vel);
hpath_length = parker_spiral_length(hd/1.496e+11,sw_vel);
figure(1)
plot(ttpk,bpath_length/1.496e+11,'b.')
hold on
plot(ttpk,hpath_length/1.496e+11,'r.')
legend('base','harmonic')
title('beam path distance')
datetick()
hold off

%dt=86400*(ttpk(1:end-1)-ttpk(2:end));
%bds=bpath_length(1:end-1)-bpath_length(2:end);
%hds=hpath_length(1:end-1)-hpath_length(2:end);
%bv=bds./dt';
%hv=hds./dt';
%figure(2)
%plot(ttpk(1:end-1),-movmean(bv,8))
%hold on
%plot(ttpk(1:end-1),-movmean(hv,8))
%legend('base','harmonic')
%title('velocity')
%datetick()
%hold off
end

