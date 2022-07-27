function ev = vel_to_ev_no(vel)

% ev = vel_to_ev(vel)
%
% Converts velocity [in m/sec] to electron kinetic energy given in eV 
% Uses relativistic formula. Electron rest energy mc^2 is subtracted.

load_plasma_constants;

%ve_nonrel = sqrt(2*energ*e_charge/m_el)
tt = m_el*c_light^2;
%ve = c_light*sqrt(1 - (tt./(energ*e_charge+tt)).^2);
ev = (1/e_charge).*((tt)./(1-(vel.^2)/(c_light^2))-tt);
ev = ev./2; % not sure why
