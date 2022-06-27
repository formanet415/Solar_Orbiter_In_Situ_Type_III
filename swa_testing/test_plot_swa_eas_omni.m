%eas1=spdfcdfread('Z:\swa\eas\2021\solo_L2_swa-eas1-nm3d-dnf_20211009T060517-20211009T120507_V01.cdf');
%info = spdfcdfinfo('Z:\swa\eas\2021\solo_L2_swa-eas1-nm3d-dnf_20211009T060517-20211009T120507_V01.cdf');
eas=spdfcdfread('Z:\swa\eas\2021\solo_L2_swa-eas2-nm3d-dnf_20211009T060517-20211009T120507_V01.cdf');
info = spdfcdfinfo('Z:\swa\eas\2021\solo_L2_swa-eas2-nm3d-dnf_20211009T060517-20211009T120507_V01.cdf');
info.Variables
info.Variables(:,1)

tdata = cdb_cdf_get_variable(eas, info, 'EPOCH');
el = cdb_cdf_get_variable(eas, info, 'SWA_EAS_ELEVATION');
el_d_u = cdb_cdf_get_variable(eas, info, 'SWA_EAS_ELEVATION_delta_upper');
el_d_l = cdb_cdf_get_variable(eas, info, 'SWA_EAS_ELEVATION_delta_lower');
az = cdb_cdf_get_variable(eas, info, 'SWA_EAS_AZIMUTH');
az_d_u = cdb_cdf_get_variable(eas, info, 'SWA_EAS_AZIMUTH_delta_upper');
az_d_l = cdb_cdf_get_variable(eas, info, 'SWA_EAS_AZIMUTH_delta_lower');
energy = cdb_cdf_get_variable(eas, info, 'SWA_EAS2_ENERGY');
data = cdb_cdf_get_variable(eas, info, 'SWA_EAS2_NM3D_DNF_Data');

for i=1:32
    for j = 1:16
        omega(i,j) = solidangle([az(i)+az_d_u(i)-180;az(i)-az_d_l(i)-180],[el(j)+el_d_u(j);el(j)-el_d_l(j)]);
        for k = 1:64
            omni(k,:) = data(k,i,j,:)*omega(i,j);
        end
    end
end

imagesc(tdata,energy(1,:),log(omni))
set(gca, 'YDir', 'normal')
datetick('keeplimits')
set(gca, 'YScale', 'log')
datetick('dd.mm HH:MM','keeplimits')
datetick('x','dd.mm HH:MM','keeplimits')
set(gca, 'YScale', 'log')
colorbar()