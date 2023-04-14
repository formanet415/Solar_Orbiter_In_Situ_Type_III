t3 = load_events();
load("epd_energies_V02.mat")
polarr = [];
out.r = [];
rall = [];
out.epd_vis=[];
out.freq=[];

time=[];
out.legacyfreq = nanmean(t3.tswf_frequency,2); %#ok<NANMEAN>
for i = 1:size(t3.year)
    year = t3.year(i);
    month = t3.month(i);
    day = t3.day(i);
    h = t3.radio_hour(i);
    m = t3.radio_minute(i);
    epd_nxt = t3.epd_nxtday(i);
    epd_h = t3.epd_hour(i);
    epd_m = t3.epd_minute(i);
    tswf_nxt = t3.tswf_nxtday(i);
    tswf_idx = t3.tswf_indexes(i,:);
    tswf_fq = t3.tswf_frequency(i,:);
    lang_nxt = t3.langmuir_nxtday(i);
    lang_h = t3.langmuir_hour(i);
    lang_m = t3.langmuir_minute(i);
    time(i) = datenum(year,month,day,h,m,0);
    [~,dt]=caadb_get_solo_orbit(time(i),4*3600); 
    out.r(i,1) = vecnorm(mean(dt,2))/1.495978707e8;
    out.epd_vis(i,1) = sum(sum(epd_energies(i,:,:)))>0;
    [ep, stat, extras] = caadb_get_solo_tds_stat(time(i),4*3600);
    lwfreq = stat.wa_med_freq(stat.wa_med_freq>10000); lwfreq = lwfreq(lwfreq<80000);
    scatter(ep,stat.wa_med_freq);
    hold on
    plot(ep, ep*0+80000)
    plot(ep, ep*0+10000)
    datetick
    hold off
    out.freq(i,1) = median(lwfreq)/1e3;
end
out.date = string(datestr(time,'yyyy/mm/dd-HH:MM'));
myTable = struct2table(out);

% save the table to a CSV file
writetable(myTable, 'events_summary.csv');

