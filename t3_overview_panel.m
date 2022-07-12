function t3_overview_panel(year, month, day, h, m, epd_nxt, epd_h, epd_m, tswf_nxt, tswf_idx, tswf_fq, lang_nxt, lang_h, lang_m)
%T3_OVERVIEW_PANEL Plotter function for individual events which analyses
%them and saves panels with the data

caa_data_paths;

rswf = tdscdf_load_l2_surv_rswf(datetime(year,month,day));
rtime0 = datenum(year, month, day, h, m, 0) - 1/48;
rtime1 = rtime0 + 2/12 - 1/48;
[~, r0] = min(abs(rswf.epoch - rtime0));
[~, r1] = min(abs(rswf.epoch - rtime1));

ridxs = r0:r1;
for i = 1:length(ridxs)
    
    srfuu = convert_to_SRF(rswf,ridxs(i));  % This takes samps_per_ch into account
    
    % commented an example of conversion to B paralell and B perpendicular
    %srfuu(1:2,:) = rswf.data(1:2,:,i);
    %[ep, b_vec, b_range, time_res, qf]=caadb_get_solo_mag(datenum(year,month,day),1,'srf');
    %b = b_vec(:,1);
    %bp = b(2:3);
    %bpn = bp/sqrt(bp'*bp);
    %bper = (srfuu'*bpn)';
    %bort = ([bpn(2),-bpn(1)]*srfuu);
    %[spz, fq, nav] = make_spectrum(srfuu(2,:), length(srfuu(2,:)), 1/rswf.samp_rate(ridxs(i)));
    
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

nsplts = 5; %number of subplots




subplot(nsplts,1,1);    % RSWF spectrogram from L2 data
ylim auto
opts.colorax = [-16,-9];
hcbar = plot_spectrogram(sps',fq,tt,opts);
xlim([rtime0,rtime1])
title(sprintf('TDS RSWF spectrogram %s',datestr(rtime0,'yyyy-mm-dd HH:MM:SS.FFF')))
cb=colorbar;
cb.Position=[0.91,0.77,0.01,0.16];
ylim manual
line([datenum(year,month,day+epd_nxt,epd_h,epd_m,0) datenum(year,month,day+epd_nxt,epd_h,epd_m,0)], get(gca, 'ylim'), 'Color', 'black');


subplot(nsplts,1,2);    % EPD STEP panel
ylim auto
if rtime0<datenum(2021,10,22)
solo_panel_epd_step_rates_spectrum(rtime0,4*60*60)
elseif rtime0>datenum(2021,10,22)
    solo_panel_epd_step_main_spectrum(rtime0,4*60*60)
end
xlim([rtime0,rtime1])

[pastt,pasden] = caadb_get_solo_swa_pas_moments(rtime0,4*60*60);
ylim manual
line([datenum(year,month,day+epd_nxt,epd_h,epd_m,0) datenum(year,month,day+epd_nxt,epd_h,epd_m,0)], get(gca, 'ylim'), 'Color', 'black');


subplot(nsplts,1,3)     % Plasma density
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
line([datenum(year,month,day+epd_nxt,epd_h,epd_m,0) datenum(year,month,day+epd_nxt,epd_h,epd_m,0)], get(gca, 'ylim'), 'Color', 'black');
%xline(datenum(year,month,day+epd_nxt,epd_h,epd_m,0));

subplot(nsplts,1,4)     % Magnetic field cone angle (tbd more about B field)
ylim auto
[ep, b_vec, b_range, time_res, qf]=caadb_get_solo_mag(rtime0,60*60*4,'rtn');
if ~isempty(ep)
    yyaxis right
    ntb=b_vec(2:3,:);
    CosTheta = max(min(ntb(1,:)/(norm(ntb)),1),-1);
    ThetaInDegrees = real(acosd(CosTheta));
    plot(ep,ThetaInDegrees,'Displayname','RTN clock angle')
    
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
    
    legend('AutoUpdate','off')
    line([datenum(year,month,day+epd_nxt,epd_h,epd_m,0) datenum(year,month,day+epd_nxt,epd_h,epd_m,0)], [0 180], 'Color', 'black');
else
    title('no MAG data')
    xlim([rtime0,rtime1])
    ylim manual
    line([datenum(year,month,day+epd_nxt,epd_h,epd_m,0) datenum(year,month,day+epd_nxt,epd_h,epd_m,0)], get(gca, 'ylim'), 'Color', 'black');
end
datetick('Keeplimits');
hold off

subplot(nsplts,1,5)     % Statistics from TSWF 
ylim auto
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

    plot(tswftt,tswffq*1e-3,'o')
    xlim([rtime0,rtime1])
    datetick('Keeplimits');
    ylim([0,100])
    ylabel('frequency [kHz]')
    title('TSWF Langmuir wave frequencies')
else
    title('No Langmuir waves recorded in TSWF')
    xlim([rtime0,rtime1])
    datetick('Keeplimits');
end
ylim manual
line([datenum(year,month,day+epd_nxt,epd_h,epd_m,0) datenum(year,month,day+epd_nxt,epd_h,epd_m,0)], get(gca, 'ylim'), 'Color', 'black');


f = gcf;
f.Position = [100 100 1700 1300];
saveas(f, ['overview plots' filesep sprintf('TYPE_III_overview_panel_%s.png',datestr(rtime0,'yyyymmdd_HHMMSS'))])
%exportgraphics(f, ['overview plots' filesep sprintf('TYPE_III_overview_panel_%s.png',datestr(rtime0,'yyyymmdd_HHMMSS'))],'Resolution','300')

end

