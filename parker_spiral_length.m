function l = parker_spiral_length(r,vsw)
%PARKER_SPIRAL_LENGTH Uses the formula by Robinson and Carins
%   Recreating the method used by Malaspina to estimate speeds of particle
%   beams.
if mean(r)<20 % astronomical units to meters if needed
    r = r.*1.496e11;
end
omega = 2.8e-6;
rw = vsw*1e3/omega;
l = 1./(2*rw) .* (r.*sqrt(r.^2+rw.^2)+(rw.^2).*log((r+sqrt(r.^2+rw.^2))./(rw)));
end

