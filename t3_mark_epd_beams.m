function t3_mark_epd_beams(rtt, epdtt, index, name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~exist('name') || isempty(name)
    name = 'epd_energies.mat';
end
if isfile(name)
    load(name)
end
rtt0 = rtt - 1/24;
rtt1 = rtt0 + 4/24;

opts.show_xlabel = 0;
if rtt < datenum(2021, 10, 22)
    solo_panel_epd_step_rates_spectrum(rtt0,4*60*60, 'electrons', opts)
elseif rtt > datenum(2021, 10, 22)
    solo_panel_epd_step_main_spectrum(rtt0,4*60*60, 'electrons', opts)
end
xlim([rtt0, rtt1])
f = gcf;
f.Position = [100 100 1300 600];
vertline(epdtt, 'black');
goodinput = false;
while ~goodinput
    [x,y] = ginput;
    marks = length(x);
    if rem(marks,2)==0
        goodinput = true;
        epd_energies(index,1,1:marks) = x;
        epd_energies(index,2,1:marks) = y;
        fprintf('saving %s, updated marks on %s/n', name, datestr(rtt))
        save(name,'epd_energies')
    elseif isempty(x)
        return
    end

end

end

