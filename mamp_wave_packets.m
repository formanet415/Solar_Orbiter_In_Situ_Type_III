function mamp_wave_packets(dn,type,index,skip_loading)
%MAMP_WAVE_PACKETS Plots mamp data and triggered snapshots showing the wave
%packets.


if ~exist('skip_loading','var') || isempty(skip_loading)
    skip_loading = 0;
end
if skip_loading == 1
    dt = dn;    
elseif strcmp(type,'tswf')
    dt = tdscdf_load_l2_surv_tswf(dn);
elseif strcmp(type,'rswf')
    dt = tdscdf_load_l2_surv_rswf(dn);
elseif strcmp(type,'sbm2')
    dt = tdscdf_load_l2_surv_tswf('E:\Dokumenty\OneDrive - Univerzita Karlova\OKF\sbm2\solo_L2_rpw-tds-sbm2-tswf-e_20211009T060650-20211009T080652_V01.cdf',0);
end
samps = dt.samples_per_ch(index);
dtt = samps/(dt.samp_rate(index)*(86400));
tt = linspace(dt.epoch(index),dt.epoch(index) + dtt,samps);
%plot(tt,dt.data(1,1:samps,index))
%xlim([dt.epoch(index) dt.epoch(index) + dtt])
%tick5 = linspace(dt.epoch(index),dt.epoch(index) + dtt,4);

%datetick('x','HH:MM:SS.FFF','keeplimits')
%xticks(tick5);
%xticklabels(datestr(tick5,'HH:MM:SS.FFF'));


caa_data_paths;
[mep, mdt, extras] = caadb_get_solo_tds_mamp(dt.epoch(index)-0.2/86400,3);
mlen=300;
plot(mep(1:mlen)+6.5e-3/86400,mdt(1,1:mlen),'-k')
hold on
plot(tt,dt.data(1,1:samps,index),'b')
xlim([dt.epoch(index) dt.epoch(index) + dtt])
ylabel('Electric field intensity [V/m]')
legend('MAMP', 'Waveform','Autoupdate','off')
plot(mep(1:mlen)+6.5e-3/86400,-mdt(1,1:mlen),'-k')
a = area(mep(1:mlen)+6.5e-3/86400, mdt(1,1:mlen)');
b = area(mep(1:mlen)+6.5e-3/86400, -mdt(1,1:mlen)');
a.FaceColor = 'g';
a.FaceAlpha = 0.3;
b.FaceColor = 'g';
b.FaceAlpha = 0.3;

title(sprintf('TDS %s snapshot from %s; samp rate %ikHz', type,datestr(dt.epoch(index)),fix(dt.samp_rate(index)/1000)))
datetick('Keeplimits')
tick5 = linspace(dt.epoch(index)-0.2/86400,dt.epoch(index) + dtt+0.8/86400,6);
xticks(tick5);
xticklabels(datestr(tick5,'HH:MM:SS.FFF'));
xlim([dt.epoch(index)-dtt*0.15 dt.epoch(index) + dtt*(1.15)])
hold off
print(gcf,sprintf('plots/mamp_waveform_%s_%s_%i.png', type, datestr(dt.epoch(index),'yyyy_mm_dd'),index),'-dpng','-r300');
end

