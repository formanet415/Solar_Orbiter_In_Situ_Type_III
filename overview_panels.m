t3 = load_events();
load("epd_energies_V02.mat")
polarr = [];
r = [];
rall = [];
if~exist("opts") || ~isfield(opts, 'autoplot_epd') || isempty(opts.autoplot_epd)
    opts.autoplot_epd = 0;
end
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
    [tmp, tr] = t3_overview_panel(year, month, day, h, m, epd_nxt, epd_h, epd_m, tswf_nxt, tswf_idx, tswf_fq, lang_nxt, lang_h, lang_m, i, epd_energies, opts);
    if ~isnan(tmp)
        sr = size(tmp);
        polarr(end+1:end+sr(1),1:16) = tmp;
        r(end+1) = tr; %#ok<SAGROW> 
    end
    rall(end+1) = tr;
end
save('polarisation_array_V03.1.mat','polarr')