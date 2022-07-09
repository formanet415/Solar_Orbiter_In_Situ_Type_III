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
    % add projection onto B field or at least convert to SRF
    srfuu = convert_to_SRF(rswf,ridxs(i));  % This takes samps_per_ch into account
    %srfuu(1:2,:) = rswf.data(1:2,:,i);
    %[ep, b_vec, b_range, time_res, qf]=caadb_get_solo_mag(datenum(year,month,day),1,'srf');
    %b = b_vec(:,1);
    %bp = b(2:3);
    %bpn = bp/sqrt(bp'*bp);
    %bper = (srfuu'*bpn)';
    %bort = ([bpn(2),-bpn(1)]*srfuu);
    [sp, fq, nav] = make_spectrum(srfuu, length(srfuu(1,:))/8, 1/rswf.samp_rate(ridxs(i)),rswf.samp_rate(ridxs(i))/2);
    %[spz, fq, nav] = make_spectrum(srfuu(2,:), length(srfuu(2,:)), 1/rswf.samp_rate(ridxs(i)));
    temp = sum(squeeze(sp)');
    sps(i,1:length(temp)) = temp;
    tt(i) = rswf.epoch(ridxs(i));
end
if r1 == length(rswf.epoch)
    li=i;
    %tbd load next day and append some stuff
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
    disp('day overflow, skipping this for now')
end

subplot(4,1,1);
opts.colorax=[-16,-9];
hcbar =plot_spectrogram(sps',fq,tt,opts);
xlim([rtime0,rtime1])
title(sprintf('TDS RSWF spectrogram %s',datestr(rtime0,'yyyy-mm-dd HH:MM:SS.FFF')))
cb=colorbar;
cb.Position=[0.91,0.77,0.01,0.16];

subplot(4,1,2);
if rtime0<datenum(2021,10,22)
solo_panel_epd_step_rates_spectrum(rtime0,4*60*60)
elseif rtime0>datenum(2021,10,22)
    solo_panel_epd_step_main_spectrum(rtime0,4*60*60)
end
xlim([rtime0,rtime1])


[pastt,pasden] = caadb_get_solo_swa_pas_moments(rtime0,4*60*60);
subplot(4,1,3)
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
l = legend();


subplot(4,1,4)
[ep, b_vec, b_range, time_res, qf]=caadb_get_solo_mag(rtime0,60*60*4,'rtn');
if ~isempty(ep)
    coneang = rad2deg(acos(b_vec(1,:)./vecnorm(b_vec)));
    plot(ep, coneang,'b','Displayname','RTN cone angle')
    datetick()
    xlim([rtime0,rtime1])
    legend()
    title('MAG RTN cone angle')
else
    title('no MAG data')
end


f = gcf;
f.Position = [100 100 1700 1000];
saveas(f, ['overview plots' filesep sprintf('TYPE_III_overview_panel_%s.png',datestr(rtime0,'yyyymmdd_HHMMSS'))])
%exportgraphics(f, ['overview plots' filesep sprintf('TYPE_III_overview_panel_%s.png',datestr(rtime0,'yyyymmdd_HHMMSS'))],'Resolution','300')

end

