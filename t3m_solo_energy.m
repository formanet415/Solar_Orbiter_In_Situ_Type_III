% Compares Malasipinas method to EPD data and fitted (by hand) particle
% beam speeds

caa_data_paths
load('polarisation_array_V02.mat');
load('epd_energies_V02.mat');

r = polarr(:,7);
vsw = polarr(:,8);
lens = parker_spiral_length(r,vsw);

rtts = polarr(:,4);
for i = 1:length(rtts)
    tt0 = rtts(i);
    rtt = rtts(i) - 1/48;
    endtt = rtt + 1/6;
    len = lens(i);
    if tt0<datenum(2021,10,22)
        solo_panel_epd_step_rates_spectrum(rtt,4*60*60,'electrons')
    elseif tt0>datenum(2021,10,22)
        solo_panel_epd_step_main_spectrum(rtt,4*60*60,'electrons')
    end
    xlim([rtt,endtt])
    f = gcf;
    f.Position = [570 469 1091 451];
    t = rtt:1/86400:endtt;
    own = t3_beam_energy(t, epd_energies);
    hold on
    plot(t, own, 'b.')
    t = tt0+1000/86400:1/86400:endtt;
    vel = t3m_beam_speed(len, tt0, t);
    men = vel_to_ev(vel);
    plot(t, men, 'r.');
    hold off
end