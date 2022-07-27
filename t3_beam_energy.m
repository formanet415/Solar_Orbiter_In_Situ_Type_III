function enb = t3_beam_energy(tt, epd_energies)
enb = 0;
load('t3_in_situ_events_V01.mat')
[~,index] = min(abs(tt-events.rtt));clear events;
ts = epd_energies(index,1,:);
ts = ts(ts~=0);
t0s = ts(1:2:end);
t1s = ts(2:2:end);

es = log(epd_energies(index,2,:));
es = es(es~=0);
e0s = es(1:2:end);
e1s = es(2:2:end);

enb = zeros(size(tt));

for i = 1:length(t0s)
    a = tt>t0s(i);
    b = tt<t1s(i);
    these = and(a,b);
    enb(these == 1) = exp(e0s(i) + (e1s(i)-e0s(i))/(t1s(i) - t0s(i))*(tt(these == 1) - t0s(i)));
end
%en0 = log(47);
%t1 = datenum(2021,10,9,7,30,0);
%en1 = log(12);
%enb = exp(en0 + (en1-en0)/(t1 - t0)*(tt - t0));
end
