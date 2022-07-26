function [polarray, r] = t3_overview_panel(year, month, day, h, m, epd_nxt, epd_h, epd_m, tswf_nxt, tswf_idx, tswf_fq, lang_nxt, lang_h, lang_m, index, epd_energies)
%T3_OVERVIEW_PANEL Plotter function for individual events which analyses
%them and saves panels with the data

caa_data_paths;
subplot = @(m,n,p) subtightplot (m,n,p,[0.04 0.04], [0.03 0.02], [0.04 0.04]);


rswf = tdscdf_load_l2_surv_rswf(datetime(year,month,day));
rtime0 = datenum(year, month, day, h, m, 0) - 1/48;
rtime1 = rtime0 + 2/12 - 1/48;
epd_time = datenum(year,month,day+epd_nxt,epd_h,epd_m,0);
lang_time = datenum(year,month,day+lang_nxt,lang_h,lang_m,0);
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

nsplts = 6;             % number of subplots
[ep, b_vec, b_range, time_res, qf]=caadb_get_solo_mag(rtime0,60*60*4,'rtn');
if isempty(b_vec)
    nsplts = 4;
end
osplts = [2 3 4 5 6 1]; % order of subplots


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
title(sprintf('TDS RSWF spectrogram %s',datestr(rtime0,'yyyy-mm-dd HH:MM:SS.FFF')))
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
temp = size(epd_energies);
hold on
for i = 1:(temp(3)/2)
    plot([epd_energies(index,1,i*2-1) epd_energies(index,1,i*2)],[epd_energies(index,2,i*2-1) epd_energies(index,2,i*2)],'b','LineWidth',3)
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

title('Plasma density measured by SolO instruments')
datetick()
ylabel('Density [cm^-3]')
xlim([rtime0,rtime1])
l = legend('AutoUpdate','off');
%ylim manual
ylim(get(gca, 'ylim'))
vertline(epd_time,'black');

[ep, b_vec, b_range, time_res, qf]=caadb_get_solo_mag(rtime0,60*60*4,'rtn');
if ~isempty(ep)
    subplot(nsplts,1,osplts(4))     % Magnetic field cone angle (tbd more about B field)
    hold off
    ylim auto

    yyaxis right
    ntb=b_vec(2:3,:);
    CosTheta = ntb(1,:)./(vecnorm(ntb));
    ThetaInDegrees = real(acosd(CosTheta)).*sign(ntb(2,:));
    color = [0.8500 0.3250 0.0980];
    plot(ep,ThetaInDegrees,'.','Color',color,'Markersize',1,'Displayname','RTN clock angle')
    ylim([-180 180])
    set(gca, 'Ycolor', color)
    ylabel('clock angle [deg]')
    %legend('Location','northeast')

    yyaxis left
    ylim manual
    ylim([0 180])
    hold on
    coneang = rad2deg(acos(b_vec(1,:)./vecnorm(b_vec)));
    plot(ep, coneang,'b-','Displayname','RTN cone angle')
    datetick('Keeplimits');
    xlim([rtime0,rtime1])
    ylabel('cone angle [deg]')
    title('MAG RTN angles')

    %if strcmp(version,'9.11.0.1769968 (R2021b)')
    %    legend('AutoUpdate','off','Location','northeast')
    %else
    %    legend('AutoUpdate','off','Location','northwest')
    %end

    vertline(epd_time,'black');
    datetick('Keeplimits');





    subplot(nsplts,1,osplts(5))     % E perpendicular / E total

    % commented an example of conversion to B paralell and B perpendicular
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
                continue
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
        ylabel('E^2_{\perp}/(E^2_{||} + E^2_{\perp})');
        vertline(epd_time,'black');
    end
end

graph = gcf;
graph.Position = [100 100 1700 1300];
saveas(graph, ['overview plots' filesep sprintf('TYPE_III_overview_panel_%s.png',datestr(rtime0,'yyyymmdd_HHMMSS'))])
%exportgraphics(f, ['overview plots' filesep sprintf('TYPE_III_overview_panel_%s.png',datestr(rtime0,'yyyymmdd_HHMMSS'))],'Resolution','300')
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
    polarray=[];
    for i = 1:length(f)
        [~, tmp] = min(abs(tswf.epoch-fep(i)));
        srfwf(1:2,:) = convert_to_SRF(tswf,tmp);
        Erms = sqrt(std(bort)^2+std(bper)^2);
        beamenerg = beam_energy(tswf.epoch(tmp),epd_energies,index);
        if beamenerg~=0
            polarray(end+1,1:8) = [f(i), Erms, beamenerg, rtime0, lang_time, epd_time, r, sw_vel];
        else
            polarr = nan;
        end
    end

else
    polarray = nan;
end



end


function enb = beam_energy(tt, epd_energies, index)
enb = 0;
ts = epd_energies(index,1,:);
ts = ts(ts~=0);
t0s = ts(1:2:end);
t1s = ts(2:2:end);

es = log(epd_energies(index,2,:));
es = es(es~=0);
e0s = es(1:2:end);
e1s = es(2:2:end);

enb = zeros(size(tt));

for i = 1:length(t0s)
    a = tt>t0s(i);
    b = tt<t1s(i);
    these = and(a,b);
    enb(these == 1) = exp(e0s(i) + (e1s(i)-e0s(i))/(t1s(i) - t0s(i))*(tt(these == 1) - t0s(i)));
end
%en0 = log(47);
%t1 = datenum(2021,10,9,7,30,0);
%en1 = log(12);
%enb = exp(en0 + (en1-en0)/(t1 - t0)*(tt - t0));
end
