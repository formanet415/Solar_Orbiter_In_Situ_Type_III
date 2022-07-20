function [beamtts,beampks] = t3_fit_electron_velocity(date, epdtt, makeplot)
%T3_FIT_ELECTRON_VELOCITY Summary of this function goes here
%   Detailed explanation goes here
date=date-(1/48);

if ~exist("makeplot") || isempty(makeplot)
    makeplot = 1;
end

if date<datenum(2021,10,22)
    [ep,el,mag,energ,extras] = caadb_get_solo_epd_step_rates(date,4*60*60);
elseif date>datenum(2021,10,22)
    [ep, ~, el, ~, mag, energ] = caadb_get_solo_epd_step_main(date,4*60*60);
end

averaged = 20;
dt=1./(el-mag);

r=size(el);
peaks=[];
index=[];
for i = 1:r(2)-averaged
    tmp = dt(:,i:i+averaged);
    ap = abs(mean(tmp,2));
    [~,loc]=findpeaks(log(ap));
    pks = energ(loc);
    for j = pks'
        index(end+1)=i;
        peaks(end+1)=j;
    end
end
if makeplot == 1
    scatter(ep(index),peaks);
    xlim([date date+4/24]);
    set(gca,'YScale','log')
    datetick('Keeplimits');
    vertline(epdtt,'black');
end



[~,idx] = min(abs(ep-epdtt));
start = mean(dt(:,idx-25:idx+25),2);
[~, loc] = findpeaks(log(abs(start)));
start = fix(median(loc));
if isnan(start)
    beampks = [];
    beamtts = [];
    return
end
beampks = [energ(start)];
beamtts = [epdtt];
cur=start;
for i=idx:r(2)-averaged
    tmp = dt(:,i:i+averaged);
    ap = abs(mean(tmp,2));
    [~,loc] = findpeaks(log(ap));
    d = loc(loc<=cur);

    if ~isempty(d)
        cur = max(d);
        beampks(end+1) = energ(cur);
        beamtts(end+1) = ep(i);
    end
end

if makeplot == 1
    hold on
    plot(beamtts,beampks,'r','LineWidth',8)

    f = gcf;
    f.Position = [100 100 900 600];

    if ~exist(['plots' filesep 'fitting'], 'dir')
        mkdir(['plots' filesep 'fitting'])
    end
    saveas(gcf, ['plots' filesep 'fitting' filesep sprintf('EPD_fit_testing_%s_50avg.png',datestr(date,'yyyymmdd_HHMMSS'))]);
    hold off
end

end

