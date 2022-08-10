t3 = load_events();
load("epd_energies_V02.mat")


version = 2;


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
    %t3_overview_panel(year, month, day, h, m, epd_nxt, epd_h, epd_m, tswf_nxt, tswf_idx, tswf_fq, lang_nxt, lang_h, lang_m, i, epd_energies)
    events.rtt(i) = datenum(year, month, day, h, m, 0);
    events.epdtt(i) = datenum(year, month, day + epd_nxt, h, m, 0);
    events.tswf_nxt(i) = tswf_nxt;
    events.tswf_idx(i,1:length(tswf_idx)) = tswf_idx;
    events.langtt(i) = datenum(year, month, day + lang_nxt, h, m, 0);
    tmp = size(epd_energies);
    events.epd_energies(i, :, 1:tmp(3)) = epd_energies(i, :, :);
end
events.epd_energies(events.epd_energies == 0) = nan;
save(sprintf('t3_in_situ_events_V%02i.mat',version),'events')