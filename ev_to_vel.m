function ve = ev_to_vel(energ)

% ve = ev_to_vel(energ)
%
% Converts electron kinetic energy given in eV to particle velocity [in m/sec]
% Uses relativistic formula. Electron rest energy mc^2 is subtracted.

load_plasma_constants;

%ve_nonrel = sqrt(2*energ*e_charge/m_el)
tt = m_el*c_light^2;
ve = c_light*sqrt(1 - (tt./(energ*e_charge+tt)).^2);

