clear; clc; close all;
%% 1
N = 128;
S = zeros(1, N);
S(8) = 2;
S(24) = 2;
S(40) = 3;
Signal = ifft(S, N);
%SignalSpec = fft(real(Signal), N);
%figure
%plot(linspace(0, 1, N), SignalSpec)

%SignalSpec = fft(imag(Signal), N);
%figure
%plot(linspace(0, 1, N), imag(SignalSpec))

%xlabel('Time')
%ylabel('Amp')
%title('Signal')
%% 2
SNR = 25;
NoisedSignal = NoiseGenerator(SNR, Signal);
figure
plot(linspace(0, 1, N), abs(NoisedSignal))
xlabel('Time')
ylabel('Amp')
title('NoisedSignal')
%% 3
Noise = NoisedSignal - Signal;
P_Signal = PowerSignal(Signal);
P_Noise = PowerSignal(Noise);
P_NoisedSignal = PowerSignal(NoisedSignal);
%% 4
SignalSpec = fft(Signal, N);
NoiseSpec = fft(Noise, N);
NoisedSignalSpec = fft(NoisedSignal, N);
P_SignalSpec = PowerSignal(SignalSpec)/N;
P_NoiseSpec = PowerSignal(NoiseSpec)/N;
P_NoisedSignalSpec = PowerSignal(NoisedSignalSpec)/N;
if ((P_SignalSpec/P_Signal >= (1 - 0.001)) && (P_NoisedSignalSpec/P_NoisedSignal >= (1 - 0.001)) && (P_NoiseSpec/P_Noise >= (1 - 0.001)))
    disp('True')
else
    disp('False')
end
%% 5
FilteredNoisedSignal = FilterSignal(NoisedSignal);
FNS = fft(FilteredNoisedSignal, N);
figure
plot(linspace(0, 1, N), abs(FNS))
figure
plot(linspace(0, 1, N), abs(FilteredNoisedSignal))
xlabel('Time')
ylabel('Amp')
title('FilteredNoisedSignal')
%% 6
SNR = 15;
SNR_NoisedSignal = zeros(1, 1e3);
SNR_FilteredNoisedSignal = zeros(1, 1e3);
for k = 1:1e5
    NoisedSignal = NoiseGenerator(SNR, Signal);
    FilteredNoisedSignal = FilterSignal(NoisedSignal);
    SNR_NoisedSignal(k) = 20 * log10(rms(Signal)/rms(NoisedSignal - Signal));
    SNR_FilteredNoisedSignal(k) = 20 * log10(rms(Signal)/rms(FilteredNoisedSignal - Signal));
end
disp(mean(SNR_FilteredNoisedSignal - SNR_NoisedSignal));
%% 7
S = [1 1 1 1 1 1];
SS = fft(S, 1e3);
figure
plot(linspace(0, 1, 1e3), real(SS))
figure
plot(linspace(0, 1, 1e3), abs(SS))