function [H_LS] = LS_CE(Y, Xp, pilot_loc, Nfft, Nps, int_opt)
% LS channel estimation function
% Inputs:
% Y = Frequency-domain received signal
% Xp = Pilot signal
% pilot_loc = Pilot location
% N = FFT size
% Nps = Pilot spacing
% int_opt = ’linear’ or ’spline’
% output:
% H_LS = LS Channel estimate
% Np = Nfft / Nps;
Np = numel(pilot_loc);
k = 1 : Np;
LS_est(k) = Y(k) ./ Xp(k); % LS channel estimation
if lower(int_opt(1)) == 'l', method = 'linear'; else method = 'spline'; end
% Linear/Spline interpolation
H_LS = interp1(pilot_loc, LS_est, linspace(1, Nfft, Nfft), method);