
ep0 = datenum(2021, 10, 9, 6, 45,0); tlen = 3600*2;

show_xlabel = 1;
show_title = 1;
time_datenum = 1;
opts_sp.colorax = 0;
opts_sp.log_y = 1;
threshold = 1;

if ~exist('species','var') || isempty(species)
    species = 'electrons';
end

if exist('opts','var') && ~isempty(opts)
    if isfield(opts,'colorax')        
		opts_sp.colorax = opts.colorax;    
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

figure(1)
%plot_spectrogram(spn, energ*1e3, ep, opts_sp);
%figure(2)
%imagesc(peaks')
%set(gca,'YDir','normal')
%figure(3)
%scatter(ttpk,dtpk)
%datetick()
try
    f1 = fit(ttpk',log(dtpk)','poly1');
catch ME
    return
end
%hold on
%plot(ttpk,exp(f1(ttpk)))
hold off
x=dtpk'-exp(f1(ttpk));
filter = abs(x)-std(x)<0;
if sum(filter)<10
    return
end
f2 = fit(ttpk(filter)',log(dtpk(filter))','poly1','Weights',exp(linspace(5,1,sum(filter))));

if ep0<datenum(2021,10,22)
    solo_panel_epd_step_rates_spectrum(ep0,2*60*60,'electrons')
elseif ep0>datenum(2021,10,22)
    solo_panel_epd_step_main_spectrum(ep0,2*60*60,'electrons')
end
ylim([1,70])
xlim([ep0,ep0+1/12])
hold on
scatter(ttpk,dtpk,'black','DisplayName','EPD STEP peaks')
plot(ttpk,exp(f2(ttpk)),'b','LineWidth',3,'DisplayName','EPD STEP beam fit');

load('polarisation_array_V03.mat');
load('epd_energies_V02.mat');
r = polarr(:,7);
vsw = polarr(:,8);
lens = parker_spiral_length(r,vsw);
rtts = polarr(:,4);
ltts = polarr(:,5);
oct9index = 72; i = oct9index;
tt0 = rtts(i);
ltt = ltts(i);
rtt = rtts(i) - 1/48;
endtt = rtt + 5/12;
len = lens(i);
f = gcf;
f.Position = [170 169 1291 651];
t = ltt:1/86400:endtt;
vel = t3m_beam_speed(len, tt0, t);
men = vel_to_ev(vel);
nrmen = nrvel_to_ev(vel);
plot(t, men*1e-3, 'r','LineWidth',3,'DisplayName','Malaspina method prediction');
legend()