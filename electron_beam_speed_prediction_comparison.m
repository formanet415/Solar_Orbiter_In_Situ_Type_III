% EPD plot comparing electron beam speed prediction to EPD data

get_fit_difference = input("pause at each plot to decide if the fit difference should be calculated? (0/1): ");
caa_data_paths
load('polarisation_array_V03.mat');
load('epd_energies_V02.mat');
load('TNR_times.mat');
difference = [];

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
    ylim([3,70])
    xlim([t0,t1])
    f = gcf;
    f.Position = [570 469 1091 651];
    
    % add hand fitted beam speed (disabled for poster)
    t = t0:10/86400:t1;
    own = [];
    for j = 1:length(t)
        own(end+1) = t3_beam_energy(t(j), epd_energies);
    end
    hold on
    %plot(t, own, 'c.')

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
    
    

    opts = [];
    % automatic fitting
    [peaktimes, fitted] = t3_auto_fit_electron_vel(TNR_time(i),3600*2,'electrons',opts,1);
    if get_fit_difference==1 && ~isempty(fitted) && ~all(isnan(men))
        good_fit = input("Should this be used to compare the fits? (0/1): ");
        if good_fit==1
            tindexes = [];
            for peak = peaktimes
                [~, tindex] = min(abs(t-peak));
                tindexes = [tindexes, tindex]; %#ok<AGROW> 
            end
            difference(1,end+1:end+length(tindexes)) = exp(fitted(peaktimes))'-1e-3*men(tindexes);
            
            disp('works? i guess')
        end
    end

    legend('Malaspina et al.,2011', 'EPD STEP peaks', 'EPD STEP beam fit')

    hold off
    
    saveas(f, ['overview plots' filesep sprintf('TYPE_III_beam_speeds_TNR_time_%s.png',datestr(t0,'yyyymmdd_HHMMSS'))])
end
if get_fit_difference == 1
    save('difference_autosave', 'difference')
end