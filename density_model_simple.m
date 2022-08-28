function [based,harmd] = density_model_simple(freq)
%DENSITY_MODEL_SIMPLE gets plasma parameters and radio wave frequency and
%returns expected distance from the Sun
%   Input
%   freq - radio frequency (Hz)
%
%   Output
%   based - distance from the sun corresponding to base frequency (m)
%   harmd - distance from the sun corresponding to first harmonic frequency (m)

me = 9.11e-31;
e = 1.60217663e-19;
eps0 = 8.8541878128e-12;
N = 1e-6*(2*pi.*[freq' freq'./2]/e).^2*me*eps0;
A = 2.75e6;
p = 2.38;
b = 0.95;
r = (A./N).^(1/p)+b;
based = r(:,1)*695500000;
harmd = r(:,2)*695500000;
end

