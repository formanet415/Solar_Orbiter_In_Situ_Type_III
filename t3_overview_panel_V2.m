function [polarray, r] = t3_overview_panel_V2(tnr_time, epd_time, tswf_idx, opts)
%T3_OVERVIEW_PANEL_V2 Improved plot. Plotter function for individual events which analyses
%them and saves panels with the data

[year, month, day, h, m, ~] = datevec(tnr_time);
if exist("opts") %#ok<EXIST>
    if isfield(opts, 'epd_nxt')
        epd_nxt = opts.epd_nxt;
    else
        epd_nxt = 0;
    end
    
    if isfield(opts, 'tswf_nxt')
        tswf_nxt = opts.tswf_nxt;
    else
        tswf_nxt = 0;
    end

    if isfield(opts, 'lang_nxt')
        lang_nxt = opts.lang_nxt;
    else
        lang_nxt = 0;
    end

    if isfield(opts, 'subplot_conf')
        subplot_conf = opts.subplot_conf;
    else
        subplot_conf = [1 1 1 1 1 1 1];
    end
    if isfield(opts, 'subplot_conf_osplts')
        osplts = opts.subplot_conf_osplts;
    else
        osplts = [1 2 3 4 5 6 7]; % default order of subplots BE ADVISED the last two magnetic dependent subplots should not be placed above any other subplots, that will break the subplots
        %[2 3 5 6 7 1 4]
    end

    if isfield(opts, 'tds_max_freq')
        tds_max_freq = opts.tds_max_freq; % In kHz
    end
    
    if isfield(opts, 'verbose')
        verbose = opts.verbose;
    end
    if isfield(opts, 'plot_duration')
        plot_duration = opts.plot_duration;
    else
    plot_duration = 4*3600;
    end
end

caa_data_paths;
% good configuration for < 7 subplots:
subplot = @(m,n,p) subtightplot (m,n,p,[0.04 0.04], [0.03 0.02], [0.04 0.04]);
% todo: add a tight configuration and move timetags to the bottom subplot
% only.

yver = version('-release'); yver = str2num(yver(1:4)); % probably not needed, was used because some function worked differently before R2019
rswf = tdscdf_load_l2_surv_rswf(datetime(year,month,day));
tswf=tdscdf_load_l2_surv_tswf(datenum(year,month,day+tswf_nxt));
rtime0 = datenum(year, month, day, h, m, 0) - 1/48;
rtime1 = rtime0 - 1/48 + plot_duration/86400;

if ~isempty(tswf)
    if ~isempty(tswf_idx)
        lang_time = tswf.epoch(tswf_idx(1));
    end
else
    lang_time = nan;
end

[~, r0] = min(abs(rswf.epoch - rtime0));
[~, r1] = min(abs(rswf.epoch - rtime1));

ridxs = r0:r1;
for i = 1:length(ridxs)
    srfuu = convert_to_SRF(rswf,ridxs(i));  % This takes samps_per_ch into account
    [sp, fq, nav] = make_spectrum(srfuu, length(srfuu(1,:))/8, 1/rswf.samp_rate(ridxs(i)),rswf.samp_rate(ridxs(i))/2);

    temp = sum(squeeze(sp)');
    sps(i,1:length(temp)) = temp; % handles changes in samp rate
    tt(i) = rswf.epoch(ridxs(i));
end
if r1 == length(rswf.epoch)
    li=i;
    % loading next day if nescessary
    disp('loading next day')
    rswf2 = tdscdf_load_l2_surv_rswf(datetime(year,month,day)+1);
    r2 = 1;
    [~, r3] = min(abs(rswf2.epoch - rtime1));
    ridxs2 = r2:r3;
    for i = 1:length(ridxs2)
        srfuu = convert_to_SRF(rswf2,ridxs2(i));
        [sp, fq, nav] = make_spectrum(srfuu, length(srfuu(1,:))/8, 1/rswf2.samp_rate(ridxs2(i)),rswf2.samp_rate(ridxs2(i))/2);
        sps(li+i,:) = sum(squeeze(sp)');
        tt(li+i) = rswf2.epoch(ridxs2(i));
    end
end

nsplts = sum(subplot_conf);             % number of subplots
[~, bv1, ~, ~, ~]=caadb_get_solo_mag(rtime0,60*60*4,'rtn');
[~, bv2, ~, ~]=caadb_get_solo_mag_LL02(rtime0,60*60*4,'rtn');
if isempty(bv1) && isempty(bv2)
    nsplts = 5; 
    disp('No MAG data, cannot plot MAG or wave polarization')
end
clear bv1 bv2

if subplot_conf(1)==1
opts.show_xlabel = 0;           % TNR spectrogram
subplot(nsplts,1,osplts(1))
solo_panel_tnr_spectrum(rtime0,plot_duration,4,opts);
xlim([rtime0,rtime1])
vertline(tnr_time,'black');
ylim([7 1000])
set(gca,'FontSize',12);
ylabel('frequency [kHz]','FontSize',12);
end

if subplot_conf(2)==1
subplot(nsplts,1,osplts(2));    % RSWF spectrogram from L2 data
ylim auto
opts.colorax = [-16,-9];
if exist('tds_max_freq') && ~isempty(tds_max_freq)
    opts.max_freq = tds_max_freq;
end
hcbar = plot_spectrogram(sps',fq*1e-3,tt,opts);
xlim([rtime0,rtime1])
title(sprintf('TDS RSWF spectrogram %s',datestr(rtime0,'yyyy-mm-dd HH:MM:SS.FFF')))
colorbar('off');
%cb=colorbar;
%cb.Position=[0.91,0.688,0.01,0.09];
hold on

if ~isnan(tswf_idx(1))
    tswftt = [];
    tswffq = [];
    for i = tswf_idx
        if isnan(i)
            continue
        end
        wf = tswf.data(1,1:tswf.samples_per_ch(i),i);
        [sp, fq, nav] = make_spectrum(wf, length(wf)/8, 1/tswf.samp_rate(i));
        tswftt(end+1) = tswf.epoch(i);
        [~, j] = max(sp);
        tswffq(end+1) = fq(j);
    end
    [ep,dt] = caadb_get_solo_tds_stat(rtime0,4*60*60);
    plot(ep,dt.wa_med_freq*1e-3,'^','Color','magenta')


    plot(tswftt,tswffq*1e-3,'r*')
    datetick('Keeplimits');
    ylabel('frequency [kHz]')
    %title('TSWF Langmuir wave frequencies')
    legend('STAT wave frequency', 'Frequency from TSWF','AutoUpdate','off')
else
    %title('No Langmuir waves recorded in TSWF')
    [ep,dt] = caadb_get_solo_tds_stat(rtime0,plot_duration);
    plot(ep,dt.wa_med_freq*1e-3,'^','Color','magenta')
    legend('STAT wave frequency','AutoUpdate','off')

    datetick('Keeplimits');
    
end

ylim manual
vertline(tnr_time,'black');
set(gca,'FontSize',12);
ylabel('frequency [kHz]','FontSize',12);
end

if subplot_conf(3)==1
subplot(nsplts,1,osplts(3));    % EPD STEP panel
ylim auto
tmpopts.show_xlabel = 0;
if rtime0<datenum(2021,10,22)
    solo_panel_epd_step_rates_spectrum(rtime0,plot_duration,'electrons',tmpopts)
elseif rtime0>datenum(2021,10,22)
    solo_panel_epd_step_main_spectrum(rtime0,plot_duration,'electrons',tmpopts)
end
xlim([rtime0,rtime1])

ylim manual
vertline(tnr_time,'black');

t3_auto_fit_electron_vel(epd_time-1/48,3600*2,'electrons',opts,1);


hold off
set(gca,'FontSize',12);
ylabel('energy [keV]','FontSize',12);
end




if subplot_conf(4)==1
opts.show_xlabel = 0;           % MAMP
subplot(nsplts,1,osplts(4))
[mamp_ep, mamp] = caadb_get_solo_tds_mamp(rtime0,plot_duration);
if ~isempty(mamp)
    plot(mamp_ep, mamp(1,:)*1e3)
    set(gca, 'YScale', 'log')
    
    title(sprintf('TDS MAMP on %s, CH1',datestr(rtime0,'yyyy-mm-dd HH:MM:SS.FFF')))
    
else
    title('No MAMP data')
end

xlim([rtime0,rtime1])
datetick('Keeplimits')
vertline(tnr_time,'black');
set(gca,'FontSize',12);
ylabel('MAMP (mV/m)','FontSize',12);
end











if subplot_conf(5)==1
subplot(nsplts,1,osplts(5))     % Plasma density
[pastt, pasden, vel_rtn] = caadb_get_solo_swa_pas_moments(rtime0,plot_duration);
ylim auto
hold off

plot(pastt,pasden,'r','DisplayName','PAS density')
if ~isempty(pastt)
    hold on
end

[biatt,biaden] = caadb_get_solo_rpw_bia_density(rtime0,plot_duration);
plot(biatt,biaden,'g','DisplayName','BIAS density')
if ~isempty(biatt)
    hold on
end

[tnrtt,tnrden]=caadb_get_solo_rpw_tnr_density(rtime0,plot_duration);
plot(tnrtt,tnrden,'b','DisplayName','TNR density')

title('Plasma density measured by SolO instruments')
datetick()
set(gca,'FontSize',12);
ylabel('Density [cm^-3]','FontSize',12);

xlim([rtime0,rtime1])
l = legend('AutoUpdate','off');
%ylim manual
ylim(get(gca, 'ylim'))
vertline(epd_time,'black');
end

[brtnep, b_vec, b_range, time_res, qf]=caadb_get_solo_mag(rtime0,plot_duration,'rtn');
if isempty(brtnep) || sum(sum(isnan(b_vec)))>length(b_vec)
    fprintf('Mag data missing, looking for low latency data\n')
    [brtnep, b_vec, b_range, qf] = caadb_get_solo_mag_LL02(rtime0,plot_duration,'rtn');

end



if ~isempty(brtnep)
    if subplot_conf(6)==1
    subplot(nsplts,1,osplts(6))     % Magnetic field cone angle
    hold off
    ylim auto

    yyaxis right
    ntb = b_vec(2:3,:);
    CosTheta = ntb(1,:)./(vecnorm(ntb));
    clockang = real(acosd(CosTheta)).*sign(ntb(2,:));
    color = [0.8500 0.3250 0.0980];
    plot(brtnep, clockang, '.', 'Color', color, 'Markersize', 1, 'Displayname', 'RTN clock angle')
    ylim([-180 180])
    set(gca, 'Ycolor', color)
    set(gca,'FontSize',12);
    ylabel('clock angle [deg]','FontSize',12);
    
    %legend('Location','northeast')

    yyaxis left
    ylim manual
    ylim([0 180])
    hold on
    coneang = rad2deg(acos(b_vec(1,:)./vecnorm(b_vec)));
    plot(brtnep, coneang,'b-','Displayname','RTN cone angle')
    datetick('Keeplimits');
    xlim([rtime0,rtime1])
    ylabel('cone angle [deg]')
    title('MAG RTN angles')
    
    %legend('AutoUpdate','off')
    vertline(tnr_time,'black');
    datetick('Keeplimits');
    end


    if subplot_conf(7)==1
    subplot(nsplts,1,osplts(7))     % E perpendicular / E total

    % conversion to B paralell and B perpendicular below
    
    if ~isnan(tswf_idx(1)) && ~isempty(b_vec)
        f = [];
        fep=[];
        for i = tswf_idx
            if isnan(i)
                continue
            end
            srfwf(1:2,:) = convert_to_SRF(tswf,i);
            [ep, srf_b_vec, b_range, time_res, qf]=caadb_get_solo_mag(tswf.epoch(i)-0.25/(86400),0.5,'srf');
            if isempty(srf_b_vec)
                fprintf('Mag data missing, looking for low latency data\n')
                [ep, srf_b_vec, b_range, qf] = caadb_get_solo_mag_LL02(rtime0,0.5,'srf');
                if isempty(srf_b_vec)
                    continue
                end
            end
            b = mean(srf_b_vec,2);
            bp = b(2:3);
            bpn = bp/sqrt(bp'*bp);
            bper = (srfwf'*bpn)'/cos(atan(b(1)/sqrt(b(2)^2+b(3)^2)));
            bort = ([bpn(2),-bpn(1)]*srfwf);
            f(end+1) = std(bort)^2/(std(bort)^2+std(bper)^2);
            fep(end+1) = tswf.epoch(i);


        end
        plot(fep, f,'b.')
        xlim([rtime0,rtime1])
        datetick('Keeplimits');
        ylim([0,1])
        title('Wave polarization  F=E^2_{\perp}/(E^2_{||} + E^2_{\perp}) [F=1 => transverse, F=0 => linear]');
        set(gca,'FontSize',12);
        ylabel('E^2_{\perp}/(E^2_{||} + E^2_{\perp})','FontSize',12);
        vertline(tnr_time,'black');
    end
    end
end

graph = gcf;
graph.Position = [100 100 1700 1450];
if opts.autoplot_epd == 1
    saveas(graph, ['overview plots' filesep sprintf('TYPE_III_overview_panel_autofit_%s.png',datestr(rtime0,'yyyymmdd_HHMMSS'))])
else
    saveas(graph, ['overview plots' filesep sprintf('TYPE_III_overview_panel_%s.png',datestr(rtime0,'yyyymmdd_HHMMSS'))])
end
close(graph)


% save polarisation, energy rms, beam speed (tbd wavenumber)

[~, pos, ~] = caadb_get_solo_orbit(rtime0, 3600);
if ~isempty(pos)
    r = pos(1)/1.496e8;
else
    r = nan;
end

tmp = vecnorm(vel_rtn);
sw_vel = mean(tmp);

if exist('f') && ~isempty(f)
    disp('_____________________________________')
    disp('THIS DAY IS BEING SAVED TO POLARR: ')
    disp(datestr(rtime0))
    disp('_____________________________________')
    polarray=[];
    for i = 1:length(f)
        [~, tswfindex] = min(abs(tswf.epoch-fep(i)));
        srfwf(1:2,:) = convert_to_SRF(tswf,tswfindex);
        aErms = sqrt(std(bort)^2+std(bper)^2);
        abeamenerg = t3_beam_energy(tswf.epoch(tswfindex),epd_energies);
        if abeamenerg == 0
            polarr = nan;
            continue
        end
        % Magnetic field statistics
        [~, magindex] = min(abs(brtnep-fep(i)));
        aclockang = mean(clockang(magindex-1:magindex+1));
        aconeang = mean(coneang(magindex-1:magindex+1));
        amagfstrength = mean(vecnorm(b_vec(:,magindex-1:magindex+1)));
        
        % Plasma density statistics
        [~, pasindex] = min(abs(pastt-fep(i)));
        [~, biaindex] = min(abs(biatt-fep(i)));
        [~, tnrindex] = min(abs(tnrtt-fep(i)));
        apas = nan; if ~isempty(pasindex) apas = pasden(pasindex); end
        abia = nan; if ~isempty(biaindex) abia = biaden(biaindex); end
        atnr = nan; if ~isempty(tnrindex) atnr = tnrden(tnrindex); end
        % custom product with priorities: 1. TNR, 2. PAS (+4%), 3. BIAS
        if ~isnan(atnr) 
            aden = atnr; 
        elseif ~isnan(apas) 
            aden = apas*1.04; 
        else
            aden = abia;
        end
        if isempty(lang_time)
	    lang_time = nan;
    	end	    
        
        % polarisation, energy rms, beam speed, radio time, ...
        % langmuir time, beam time, distance from the Sun in AU, ...
        % solar wind velocity, clock angle, cone angle, ...
        % magnetic field strength, TNR density, PAS density, BIAS density,
        % combined density, timetag for cross-checking
        polarray(end+1,1:16) = [f(i), aErms, abeamenerg, rtime0, ...
            lang_time, epd_time, r, sw_vel, aclockang, aconeang, ...
            amagfstrength, atnr, apas, abia ...
            aden, fep(i)];
        
    end

else
    polarray = nan;
end

%close(graph)
end



