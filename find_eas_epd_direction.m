% Using this script I want to find out which bins from EAS1 and EAS2 are
% pointing in the general direction that is the same as the EPD sensor (its
% 15 pixels).

coslimit = 0.99; % corresponds to 66% simulated angular response of the individual STEP pixels (Rodrigues-Pacheco et al. 2020, Fig 11)

epd = cdf_load_tswf('Z:\epd\step\rates\2021\solo_L2_epd-step-rates_20211009_V01.cdf');
%epdvec = -epd.XYZ.data; % XYZ is the flow direction opposite to the view
%that would be if only considering the integral flux, lets be a bit more
%precise

epdvecs = -epd.XYZ_Sectors.data; % the flow direction opposite to the view

% EAS1 directions
[~, ~, elevation1, ~, azimuth1, ~, extras1] = caadb_get_solo_swa_eas1_nm3d_dnf(datenum(2021,10,9,6,30,0), 3600*4);
elevations1 = repmat(elevation1(:,1),1,32);
azimuths1 = repmat(azimuth1(:,1),1,16)';
r1 = zeros(16,32,3);
[r1(:,:,1),r1(:,:,2),r1(:,:,3)] = sph2cart(deg2rad(azimuths1), deg2rad(elevations1), 1*ones(16,32));
r1 = reshape(r1, 16*32,3);

figure;
hold on;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
% Plot arrows for each vector
vectors = extras1.eas_to_srf(:,:,1)*r1';
eas1mask=zeros(1,16*32);
for i = 1:size(vectors, 2)
quiver3(0, 0, 0, vectors(1, i), vectors(2, i), vectors(3, i), 'r', 'LineWidth', 0.5);
if any([vectors(1, i), vectors(2, i), vectors(3, i)]*epdvecs'>coslimit)
    eas1mask(i)=1;
    quiver3(0, 0, 0, vectors(1, i), vectors(2, i), vectors(3, i), 'y', 'LineWidth', 1);
end
end
eas1mask = reshape(eas1mask,16,32);
% Set the view angle and limits
view(3);
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
% Add a legend
legend('EAS1');
% Adjust the figure properties as needed (e.g., background color, title, etc.)
set(gcf, 'Color', 'white');
title('Pixels locations');




% EAS2 directions
[~, ~, elevation2, ~, azimuth2, ~, extras2] = caadb_get_solo_swa_eas2_nm3d_dnf(datenum(2021,10,9,6,30,0), 3600*4);
elevations2 = repmat(elevation2(:,1),1,32);
azimuths2 = repmat(azimuth2(:,1),1,16)';
r2 = zeros(16,32,3);
[r2(:,:,1),r2(:,:,2),r2(:,:,3)] = sph2cart(deg2rad(azimuths2), deg2rad(elevations2), 1*ones(16,32));
r2 = reshape(r2, 16*32,3);

axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
% Plot arrows for each vector
vectors = extras2.eas_to_srf(:,:,1)*r2';
eas2mask=zeros(1,16*32);
for i = 1:size(vectors, 2)
quiver3(0, 0, 0, vectors(1, i), vectors(2, i), vectors(3, i), 'g', 'LineWidth', 0.5);
if any([vectors(1, i), vectors(2, i), vectors(3, i)]*epdvecs'>coslimit)
    eas2mask(i)=1; 
    quiver3(0, 0, 0, vectors(1, i), vectors(2, i), vectors(3, i), 'y', 'LineWidth', 1);
end
end
eas2mask = reshape(eas2mask,16,32);
% Set the view angle and limits
view(3);
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
% Add a legend
legend('EAS1');
% Adjust the figure properties as needed (e.g., background color, title, etc.)
set(gcf, 'Color', 'white');
title('Pixels locations');


quiver3(0, 0, 0, epdvec(1), epdvec(2), epdvec(3), 'b', 'LineWidth', 1.5);
