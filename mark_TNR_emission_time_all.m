function mark_TNR_emission_time_all(name)
t3 = load_events();
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
    epdtt = datenum(year, month, day+epd_nxt, epd_h, epd_m, 0);
    rtime0 = datenum(year, month, day, h, m, 0);
    mark_TNR_emission_time(rtime0, i, name)
end
load(name)
t3.TNR_time = TNR_time; % earliest time when the type III radio emission was detected
t3.emission_time = emission_time; % the delay due to finite speed of light subtracted