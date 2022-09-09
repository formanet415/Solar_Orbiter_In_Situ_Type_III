t3 = load_events();
caa_data_paths
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
    epd_logged_trigger = datenum(year,month,day+epd_nxt,epd_h,epd_m,0);
    t3_auto_fit_electron_vel(datenum(year, month, day, h, m,00),3600*2);

    f = gcf;
    %f.Position = [100 100 1700 1300];
    datetick('Keeplimits')
    saveas(f, ['overview plots' filesep sprintf('fit_testing_TYPE_III_overview_panel_%s.png',datestr(datenum(year, month, day, h, m,00),'yyyymmdd_HHMMSS'))])
end