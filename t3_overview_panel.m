function [polarray, r] = t3_overview_panel(year, month, day, h, m, epd_nxt, epd_h, epd_m, tswf_nxt, tswf_idx, tswf_fq, lang_nxt, lang_h, lang_m, index, epd_energies, opts)
%T3_OVERVIEW_PANEL Plotter function for individual events which analyses
%them and saves panels with the data

caa_data_paths;
subplot = @(m,n,p) subtightplot (m,n,p,[0.024 0.024], [0.03 0.04], [0.03 0.03]);

yver = version('-release'); yver = str2num(yver(1:4));
rswf = tdscdf_load_l2_surv_rswf(datetime(year,month,day));
rtime0 = datenum(year, month, day, h, m, 0) - 1/48;
rtime1 = rtime0 + 2/12 - 1/48;
epd_time = datenum(year,month,day+epd_nxt,epd_h,epd_m,0);
if yver<2019
    if strcmp(lang_h{1},'nan')
        lang_time = nan;
    else
        lang_time = datenum(year,month,day+lang_nxt,str2num(lang_h{1}),str2num(lang_m{1}),0);
    end
else
    lang_time = datenum(year,month,day+lang_nxt,lang_h,lang_m,0);
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

nsplts = 7;             % number of subplots
[~, bv1, ~, ~, ~]=caadb_get_solo_mag(rtime0,60*60*4,'rtn');
[~, bv2, ~, ~]=caadb_get_solo_mag_LL02(rtime0,60*60*4,'rtn');
if isempty(bv1) && isempty(bv2)
    nsplts = 5;
end
clear bv1 bv2
osplts = [2 3 5 6 7 1 4]; % order of subplots


opts.show_xlabel = 0;           % MAMP
subplot(nsplts,1,osplts(7))
[mamp_ep, mamp, tmp] = caadb_get_solo_tds_mamp(rtime0,4*3600);
if ~isempty(mamp)
    plot(mamp_ep, mamp(1,:)*1e3)
    set(gca, 'YScale', 'log')
    
    title(sprintf('TDS MAMP CH%i %s', tmp.channel_cfg(1,1), datestr(rtime0,'yyyy-mm-dd HH:MM:SS.FFF')),'FontSize',12 )
    ylabel('MAMP (mV/m)')
else
    title('No MAMP data')
end

xlim([rtime0,rtime1])
datetick('Keeplimits')
vertline(epd_time,'black');


opts.show_xlabel = 0;           % TNR spectrogram
subplot(nsplts,1,osplts(6))
solo_panel_tnr_spectrum(rtime0,4*3600,4,opts);
xlim([rtime0,rtime1])
vertline(epd_time,'black');
ylim([7 1000])



subplot(nsplts,1,osplts(1));    % RSWF spectrogram from L2 data
ylim auto
opts.colorax = [-16,-9];
hcbar = plot_spectrogram(sps',fq*1e-3,tt,opts);
xlim([rtime0,rtime1])
title(sprintf('TDS RSWF spectrogram %s',datestr(rtime0,'yyyy-mm-dd HH:MM:SS.FFF')),'FontSize',12)
colorbar('off');
%cb=colorbar;
%cb.Position=[0.91,0.688,0.01,0.09];
hold on
tswf=tdscdf_load_l2_surv_tswf(datenum(year,month,day+tswf_nxt));
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
    [ep,dt] = caadb_get_solo_tds_stat(rtime0,4*60*60);
    plot(ep,dt.wa_med_freq*1e-3,'^','Color','magenta')
    legend('STAT wave frequency','AutoUpdate','off')

    datetick('Keeplimits');
end

ylim manual
vertline(epd_time,'black');



subplot(nsplts,1,osplts(2));    % EPD STEP panel
ylim auto
tmpopts.show_xlabel = 0;
if rtime0<datenum(2021,10,22)
    solo_panel_epd_step_rates_spectrum(rtime0,4*60*60,'electrons',tmpopts)
elseif rtime0>datenum(2021,10,22)
    solo_panel_epd_step_main_spectrum(rtime0,4*60*60,'electrons',tmpopts)
end
xlim([rtime0,rtime1])

ylim manual
vertline(epd_time,'black');
if opts.autoplot_epd == 0
    temp = size(epd_energies);
    hold on
    for i = 1:(temp(3)/2)
        plot([epd_energies(index,1,i*2-1) epd_energies(index,1,i*2)],[epd_energies(index,2,i*2-1) epd_energies(index,2,i*2)],'b','LineWidth',3)
    end
end
if opts.autoplot_epd == 1
    t3_auto_fit_electron_vel(epd_time-1/48,3600*2,'electrons',opts,1);
end


hold off



subplot(nsplts,1,osplts(3))     % Plasma density
[pastt, pasden, vel_rtn] = caadb_get_solo_swa_pas_moments(rtime0,4*60*60);
ylim auto
hold off

plot(pastt,pasden,'r','DisplayName','PAS density')
if ~isempty(pastt)
    hold on
end

[biatt,biaden] = caadb_get_solo_rpw_bia_density(rtime0,4*60*60);
plot(biatt,biaden,'g','DisplayName','BIAS density')
if ~isempty(biatt)
    hold on
end

[tnrtt,tnrden]=caadb_get_solo_rpw_tnr_density(rtime0,4*60*60);
plot(tnrtt,tnrden,'b','DisplayName','TNR density')

title('Plasma density measured by SolO instruments','FontSize',12)
datetick()
ylabel('Density [cm^-3]')
xlim([rtime0,rtime1])
l = legend('AutoUpdate','off');
%ylim manual
ylim(get(gca, 'ylim'))
vertline(epd_time,'black');

[brtnep, b_vec, b_range, time_res, qf]=caadb_get_solo_mag(rtime0,60*60*4,'rtn');
if isempty(brtnep) || sum(sum(isnan(b_vec)))>length(b_vec)
    fprintf('Mag data missing, looking for low latency data\n')
    [brtnep, b_vec, b_range, qf] = caadb_get_solo_mag_LL02(rtime0,60*60*4,'rtn');
end
if ~isempty(brtnep)
    subplot(nsplts,1,osplts(4))     % Magnetic field cone angle
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
    ylabel('clock angle [deg]')
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
    title('MAG RTN angles','FontSize',12)
    
    %legend('AutoUpdate','off')
    vertline(epd_time,'black');
    datetick('Keeplimits');



    subplot(nsplts,1,osplts(5))     % E perpendicular / E total

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
                [ep, srf_b_vec, b_range, qf] = caadb_get_solo_mag_LL02(rtime0,60*60*4,'srf');
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
        title('Wave polarization  F=E^2_{\perp}/(E^2_{||} + E^2_{\perp}) [F=1 => transverse, F=0 => linear]','FontSize',12);
        ylabel('E^2_{\perp}/(E^2_{||} + E^2_{\perp})');
        vertline(epd_time,'black');
    end
end

graph = gcf;
graph.Position = [100 100 1700 1450];
% Get handles of all subplots
h = get(gcf, 'children');

% Loop over subplots and remove XTickLabels except for the bottom one
for i = 1:length(h)
    if strcmp(get(h(i), 'type'), 'axes')  % check if this is a subplot
        if i ~= 1  % check if this is the bottom subplot
            set(h(i), 'XTickLabel', [])
        end
    end
end
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
            abeamenerg = nan;
        end
        % Magnetic field statistics
        [~, magindex] = min(abs(brtnep-fep(i)));
        if isnan(magindex) || magindex == 1 
            aclockang = nan;
            aconeang = nan;
            amagfstrength = nan;
        else
            aclockang = mean(clockang(magindex-1:magindex+1));
            aconeang = mean(coneang(magindex-1:magindex+1));
            amagfstrength = mean(vecnorm(b_vec(:,magindex-1:magindex+1)));
        end
        
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
        
        ratestotalenerg = nan;
        maintotalenerg = nan;

        if rtime0<datenum(2021,10,22)
            [epdep,int_flux, mag_flux, energ_cent, extras] = caadb_get_solo_epd_step_rates(tswf.epoch(tswfindex)-30/86400,60);
            if ~isempty(epdep) && 1 < length(epdep)
                denerg(2:length(energ_cent)) = energ_cent(2:end)-energ_cent(1:end-1);
                denerg(1)=denerg(2);

                dtepd = (epdep(2:end)-epdep(1:end-1))/86400;
                dtepd(end+1)=dtepd(end);
                ratestotalenerg = (dtepd*int_flux'*denerg')/sum(dtepd);
            end
        elseif rtime0>datenum(2021,10,22)
            [epdep, int_flux, int_flux_avg, mag_flux, mag_flux_avg, energ_cent, extras] = caadb_get_solo_epd_step_main(tswf.epoch(tswfindex)-30/86400,60);
            if ~isempty(epdep) && 1 < length(epdep)
                denerg(2:length(energ_cent)) = energ_cent(2:end)-energ_cent(1:end-1);
                denerg(1)=denerg(2);

                dtepd = (epdep(2:end)-epdep(1:end-1))/86400;
                dtepd(end+1)=dtepd(end);
                maintotalenerg = (dtepd*int_flux_avg'*denerg')/sum(dtepd);
            end
        end
        

        

        % polarisation, energy rms, beam speed, radio time, ...
        % langmuir time, beam time, distance from the Sun in AU, ...
        % solar wind velocity, clock angle, cone angle, ...
        % magnetic field strength, TNR density, PAS density, BIAS density,
        % combined density, timetag for cross-checking
        polarray(end+1,1:18) = [f(i), aErms, abeamenerg, rtime0, ...
            lang_time, epd_time, r, sw_vel, aclockang, aconeang, ...
            amagfstrength, atnr, apas, abia ...
            aden, fep(i), ratestotalenerg, maintotalenerg];
        
    end

else
    polarray = nan;
end

%close(graph)
end



