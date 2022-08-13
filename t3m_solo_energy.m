% Compares Malasipinas method to EPD data and fitted (by hand) particle
% beam speeds

caa_data_paths
load('polarisation_array_V03.mat');
load('epd_energies_V02.mat');

r = polarr(:,7);
vsw = polarr(:,8);
lens = parker_spiral_length(r,vsw);

rtts = polarr(:,4);
ltts = polarr(:,5);
for i = 1:length(rtts)
    tt0 = rtts(i);
    ltt = ltts(i);
    rtt = rtts(i) - 1/48;
    endtt = rtt + 1/6;
    len = lens(i);
    if tt0<datenum(2021,10,22)
        solo_panel_epd_step_rates_spectrum(rtt,4*60*60,'electrons')
    elseif tt0>datenum(2021,10,22)
        solo_panel_epd_step_main_spectrum(rtt,4*60*60,'electrons')
    end
    ylim([1,70])
    xlim([rtt,endtt])
    f = gcf;
    f.Position = [570 469 1091 651];
    t = rtt:10/86400:endtt;
    own = [];
    for j = 1:length(t)
        own(end+1) = t3_beam_energy(t(j), epd_energies);
    end
    hold on
    plot(t, own, 'b.')
    t = ltt:1/86400:endtt;
    vel = t3m_beam_speed(len, tt0, t);
    men = vel_to_ev(vel);
    nrmen = nrvel_to_ev(vel);
    plot(t, men*1e-3, 'g.');
    plot(t, nrmen*1e-3, 'y.');
    set(gca, "YScale", "log")
    
    vertline(polarr(i,5));
    vertline(polarr(i,4));
    vertline(polarr(i,6));
    legend('hand fitted', 'malaspina', 'malaspina nonrelativistic')
    hold off
    
    saveas(f, ['overview plots' filesep sprintf('TYPE_III_beam_speeds_%s.png',datestr(tt0,'yyyymmdd_HHMMSS'))])
end