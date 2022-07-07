function E = convert_to_SRF(data, index)
%CONVERT_TO_SRF Summary of this function goes here
%   Detailed explanation goes here
% always convert only into SRF Y-Z plane!
if ~exist("index", 'var')
    index = 1;
end

E = zeros(2,data.samples_per_ch(index));
tds_mode = convertCharsToStrings(char(data.tds_config_label(:,index)));
if contains(tds_mode, 'SE1')
    % Pachenko's antenna angle
    pacang = 125.;
    V1=[0, 1];
    % Pachenko monopole antennas
    V2=[sind(pacang), cosd(pacang)];
    V3=[-sind(pacang), cosd(pacang)];
    % SE1 TDS mode
    M(1,:) = V1; %CH1
    M(2,:) = V2; %CH2
else
    % ANT12 158.1 deg, ANT13 -158.2 deg from SRF-Z in Y-Z plane
    pacang = 158.1;
    ant21= [sind(pacang), cosd(pacang)]; %E-field in the same sense.
    pacang = -158.2;
    ant13= -1. * [sind(pacang), cosd(pacang)]; %ant31 then -1. to ant13
    M(1,:) = ant13; %
    M(2,:) = ant21; %
end
nsamp = double(data.samples_per_ch(index));
ww = data.data(:,1:nsamp,index);
% projection: E = MAT(ANT->SRF) * V; where MAT(2,2) and V is observed field
M = inv(M);
E = M*ww(1:2,:); % transformation into SRF (Y-Z) this is (V/m)


