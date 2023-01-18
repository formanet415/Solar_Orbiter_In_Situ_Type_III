function val = mark_TNR_emission_time(rtt, index, name)
% mark_TNR_emission_time: 
%    This function plots the TNR spectrogram and recieves input from user to
%   mark the earliest time when the type III radio emission occurs. We save
%   this timetag along with the estimated emission time (subtracting the
%   speed of light delay) to a .mat file chosen by the name parameter. 
%    If no time is selected function does nothing and returns zero.
%    If you missclick, press any letter while pointing at the graph to retry
%   selection.
%    Each selection is ended by pressing enter (using ginput).
if ~exist('name') || isempty(name)
    name = 'TNR_times_autosave.mat';
end
if isfile(name)
    load(name)
end
rtt0 = rtt - 1/24;
rtt1 = rtt0 + 4/24;

opts.show_xlabel = 0;   % Plotting the TNR spectrum of event
for i=1:3
    subtightplot(3,1,i)
    solo_panel_tnr_spectrum(rtt0, 4*60*60, i, opts);
    vertline(rtt, 'black');
end
xlim([rtt0, rtt1])
f = gcf;
f.Position = [100 100 1800 450];


% Getting sc position at the current date
[~,tmp] = caadb_get_solo_orbit(rtt,1e6);
r = tmp(1,1)*1e3; % distance in meters
time_delay = r/299792458;

goodinput = false;
while ~goodinput
    [x,~,button] = ginput;
    marks = length(x);
    if all(button==1)
        goodinput = true;
        TNR_time(index:index+(marks-1)) = x;
        emission_time(index:index+(marks-1)) = x - time_delay/(86400);
        fprintf('saving %s, updated times on %s, indexes %i-%i \n', name, datestr(rtt), index, index+marks-1)
        save(name,'TNR_time','emission_time')
        val = marks;
    elseif isempty(x)
        val = 0;
        return
    else
        disp(['Incorrect input, try again. Use the mouse to select the ' ...
            'earliest time at which the radio emission appears. If you ' ...
            'missclick press C while pointing at the graph.'])
    end

end

end

