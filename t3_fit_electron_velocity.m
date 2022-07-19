function t3_fit_electron_velocity(date, epdtt)
%T3_FIT_ELECTRON_VELOCITY Summary of this function goes here
%   Detailed explanation goes here
date=date-(1/48);
[ep,el,mag,energ,extras]=caadb_get_solo_epd_step_rates(date,4*60*60);
something = 'TBD';
averaged = 50;
dt=el-mag;

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
scatter(ep(index),peaks);


xlim([date date+4/24]);
set(gca,'YScale','log')
datetick('Keeplimits');
vertline(epdtt,'black');




[~,idx] = min(abs(ep-epdtt));
start = mean(dt(:,idx-25:idx+25),2);
[~, loc] = findpeaks(log(abs(start)));
start = fix(median(loc));
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
hold on
plot(beamtts,beampks,'r','LineWidth',8)




saveas(gcf, ['plots' filesep 'fitting' filesep sprintf('EPD_fit_testing_%s_50avg.png',datestr(date,'yyyymmdd_HHMMSS'))]);

end

