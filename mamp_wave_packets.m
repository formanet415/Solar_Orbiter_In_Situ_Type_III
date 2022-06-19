function mamp_wave_packets(dn,type,index,skip_loading)
%MAMP_WAVE_PACKETS Summary of this function goes here
%   Detailed explanation goes here


if ~exist('skip_loading','var') || isempty(skip_loading)
    skip_loading = 0;
end
if skip_loading == 1
    dt = dn;    
elseif strcmp(type,'tswf')
    dt = tdscdf_load_l2_surv_tswf(dn);
elseif strcmp(type,'rswf')
    dt = tdscdf_load_l2_surv_rswf(dn);
end
t = tiledlayout(1,1);
ax1 = axes(t);

plot(dt.data(1,1:dt.samples_per_ch(index),index))
xlim([1 dt.samples_per_ch(100)])
end
% it crashes here, broken loader function?
[mep, mdt, extras] = caadb_get_solo_tds_mamp(dt.epoch(index),1);
ax2 = axes(t);
plot(mamp(2,1:10),'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
xlim([0.5,7.5])
ylim([-0.02,0.02])
end

