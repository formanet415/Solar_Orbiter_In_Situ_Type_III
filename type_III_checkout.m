function type_III_checkout(date, tlen)
%TYPE_III_CHECKOUT Makes plots for given datetime and saves them in a
%subfolder
caa_data_paths;

if ~exist('tlen','var') || isempty(tlen)
    tlen = 60*60*3;
end
if date<datenum(2021,10,23)
    solo_panel_epd_step_rates_spectrum(date,tlen)
else
    solo_panel_epd_step_main_spectrum(date,tlen)
end
y=year(date);
m=month(date);
d=day(date);
folder = ['events' filesep sprintf('%i%i%i',y,m,d)];
if ~exist(folder, 'dir')
    mkdir(folder)
end
saveas(gcf, [folder, filesep, sprintf('EPD_STEP_%s.png',datestr(date,'yyyymmdd_HHMMSS'))])
saveas(gcf, ['events' filesep 'all', filesep, sprintf('EPD_STEP_%s.png',datestr(date,'yyyymmdd_HHMMSS'))])
end

