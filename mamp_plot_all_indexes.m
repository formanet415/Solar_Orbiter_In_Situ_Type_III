function mamp_plot_all_indexes(dn,type,skip_loading)
%MAMP_PLOT_ALL_INDEXES Summary of this function goes here
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
elseif strcmp(type,'sbm2')
    dt = tdscdf_load_l2_surv_tswf('C:\Users\tform\OneDrive - Univerzita Karlova\OKF\sbm2\solo_L2_rpw-tds-sbm2-tswf-e_20211009T060650-20211009T080652_V01.cdf',0);
end
for i = 1:length(dt.epoch)
    mamp_wave_packets(dt,type,i,1)
end

end

