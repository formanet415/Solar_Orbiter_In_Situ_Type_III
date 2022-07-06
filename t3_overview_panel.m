function t3_overview_panel(year, month, day, h, m, epd_nxt, epd_h, epd_m, tswf_nxt, tswf_idx, tswf_fq, lang_nxt, lang_h, lang_m)
%T3_OVERVIEW_PANEL Summary of this function goes here
%   Detailed explanation goes here

rswf = tdscdf_load_l2_surv_rswf(datetime(year,month,day));
rtime0 = datenum(year, month, day, h, m, 0) - 1/48;
rtime1 = rtime0 + 2/12 - 1/48;
[~, r0] = min(abs(rswf.epoch - rtime0));
[~, r1] = min(abs(rswf.epoch - rtime1));

ridxs = r0:r1;
for i = 1:length(ridxs)
    % add projection onto B field or at least convert to SRF
    srfuu = convert_to_SRF(rswf,ridxs(i));  % This takes samps_per_ch into account
    [sp, fq, nav] = make_spectrum(uu, wnd_size, dt, maxfq, ovrlp);
end
[sp, fq, nav] = make_spectrum(uu, wnd_size, dt, maxfq, ovrlp);
disp('xd')

if r1 == length(rswf.epoch)
    %tbd load next day and append some stuff
    disp('day overflow, skipping this for now')
end



end

