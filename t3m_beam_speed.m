function v = t3m_beam_speed(len, tt0, tt)
%t3m_beam_speed Finds the beam speed along the Parker spiral 
%   len - length should be obtained by parker_spiral_length
%   tt0 - time of beggining of the radio emission
%   tt  - time at which to find the current speed

time = (tt-tt0)*86400;
v = len./time;
end

