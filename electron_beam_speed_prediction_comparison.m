% EPD plot comparing electron beam speed prediction to EPD data

caa_data_paths
load('polarisation_array_V03.mat');
load('epd_energies_V02.mat');
load('TNR_times.mat');


%lens = parker_spiral_length(r,vsw);

for i = 1:length(TNR_time)
    if TNR_time(i) == 0
        continue
    end
    % assign time values 
    te = emission_time(i);
    t0 = te - 1/48;
    t1 = t0 + 1/6;

    % plot epd spectrum
    if te<datenum(2021,10,22)
        solo_panel_epd_step_rates_spectrum(t0,4*60*60,'electrons')
    elseif te>datenum(2021,10,22)
        solo_panel_epd_step_main_spectrum(t0,4*60*60,'electrons')
    end
    ylim([1,70])
    xlim([t0,t1])
    f = gcf;
    f.Position = [570 469 1091 651];
    
    % add hand fitted beam speed
    t = t0:10/86400:t1;
    own = [];
    for j = 1:length(t)
        own(end+1) = t3_beam_energy(t(j), epd_energies);
    end
    hold on
    plot(t, own, 'b.')

    % get data for prediction based on time since emission
    [~, ~, vel_rtn] = caadb_get_solo_swa_pas_moments(te,3*60*60);
    tmp = vecnorm(vel_rtn);
    vsw = mean(tmp); % solar wind velocity
    [~,tmp] = caadb_get_solo_orbit(te,1e6);
    r = tmp(1,1)*1e3; % distance in meters
    len = parker_spiral_length(r,vsw);
    
    t = TNR_time(i):1/86400:t1;
    vel = t3m_beam_speed(len, te, t);
    men = vel_to_ev(vel);
    %nrmen = nrvel_to_ev(vel);
    plot(t, men*1e-3, 'g.');
    %plot(t, nrmen*1e-3, 'y.');
    set(gca, "YScale", "log")
    
    legend('hand fitted', 'Malaspina et al.,2011')
    hold off
    
    saveas(f, ['overview plots' filesep sprintf('TYPE_III_beam_speeds_TNR_time_%s.png',datestr(t0,'yyyymmdd_HHMMSS'))])
end