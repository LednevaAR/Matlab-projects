clear; clc; close all;
%%
LH = 2^11 - 1;
%%
coeff = 4005;
m_seq = RSLOS(LH, coeff);
disp(sum(m_seq)) % число единиц
disp(abs(sum(m_seq - 1))) % число нулей
m_seq1 = mapping(m_seq, "BPSK");
AutoCorrm = AutoCorr(m_seq1);
BN = linspace(1, numel(AutoCorrm), numel(AutoCorrm));
f = figure();
plot(BN, AutoCorrm);
xlabel('Bit number');
ylabel('AutoCorrelation function');
title('AutoCorrelation function of m-seq');
saveas(f, 'ACF_m_seq.fig')
%%
Ts = 1 / 10;
SNR = 20;
Signal = mapping(m_seq, "QPSK");
Signal1 = interp(Signal, 2);
Signal1 = circshift(Signal1, 1);
Signal1 = downsample(Signal1, 2);
%Signal1 = ifft(fft(Signal) .* exp(-1j * 2 * pi * STO * (1 / Ts) * (0 : length(fft(Signal)) - 1) / length(fft(Signal))));
Signal1 = awgn(Signal1, SNR, 'measured');
ted_output1 = TED(Signal1, Ts);
e = ted_output1(end);
%h(1) = - (1 / 6) * (e - 1) * (e - 2) * (e - 3);
%h(2) = - (1 / 2) * e * (e - 2) * (e - 3);
%h(3) = - (1 / 2) * e * (e - 1) * (e - 3);
%h(4) = - (1 / 6) * e * (e - 1) * (e - 2);
h(1) = 1 - e;
h(2) = e;
Signal2 = conv(Signal1, h, 'same');
f = figure;
scatter(real(Signal1), imag(Signal1));
f = figure;
scatter(real(Signal2), imag(Signal2));
%%
Ts = 1 / 1000;
SNR = 20;
SCO = 4;
Signal = mapping(m_seq, "QPSK");
Signal1 = resample(Signal, SCO, 1);
Signal1 = awgn(Signal1, SNR, 'measured');
ted_output2 = TED(Signal1, Ts);
%h(1) = - (1 / 6) * (e - 1) * (e - 2) * (e - 3);
%h(2) = - (1 / 2) * e * (e - 2) * (e - 3);
%h(3) = - (1 / 2) * e * (e - 1) * (e - 3);
%h(4) = - (1 / 6) * e * (e - 1) * (e - 2);
e = ted_output2(end);
h(1) = 1 - e;
h(2) = e;
Signal2 = conv(Signal1, h, 'same');
f = figure;
scatter(real(Signal1), imag(Signal1));
f = figure;
scatter(real(Signal2), imag(Signal2));
%%
%По поводу символьной синхронизации несколько замечаний.
%Моделирование SCO проще производить так, как я указал в замечаниях в канале курса.
%По поводу твоего TED: пример на сайте wirelesspi дан для действительного сигнала,
%а ты пытаешься применить его к комплексному. Что такое sign от комплексного числа -
%- не совсем понятно. 
%Можно для простоты ограничиться BPSK, если собираешься работать с этим TED.
%Zero-crossing предназначен для работы в петле обратной связи, 
%он выдаёт просто что-то пропорциональное (с учётом знака) временному смещению.
%Чтобы всё заработало, надо доделать петлю.
