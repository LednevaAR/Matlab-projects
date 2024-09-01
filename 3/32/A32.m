clear; clc; close all;
%%
Fs = 44000;
T = 3;
A = 3;
f = 1000;
Time = 0 : 1/Fs : T;
audio = A * sin(Time * 2 * pi * f);
figure
plot(Time, audio)
xlabel('Time')
ylabel('Amp')
title('Signal')
ylim([-A A])
f = (0 : length(audio) - 1) * Fs / length(audio);
i = f >= 0 & f <= Fs/2;
f1 = f(i);
Spec = fft(audio);
figure
plot(f1, abs(Spec(i)))
xlabel('Frequency')
ylabel('Amp')
title('Spec')
%audiowrite('task2.wav', audio/A, Fs)
%%
C = 2;
audio(audio > C) = C;
audio(audio < -C) = -C;
figure
plot(Time, audio)
xlabel('Time')
ylabel('Amp')
title('ClippedSignal')
ylim([-A A])
f = (0 : length(audio) - 1) * Fs / length(audio);
ff = f - Fs/2;
i = f >= 0 & f <= Fs/2;
j = f > Fs/2;
f1 = f(i);
Spec = fft(audio);
f = (0 : length(Spec) - 1) * Fs / (length(audio));
ff = f - Fs/2;
i = f >= 0 & f <= Fs/2;
j = f > Fs/2;
Spec1 = cat(2, Spec(j), Spec(i));
figure
plot(ff, abs(Spec1))
xlabel('Frequency')
ylabel('Amp')
title('ClippedSpec')
xlim([-Fs/2 Fs/2])
ylim([0 1000])
%audiowrite('task2new.wav', audio/A, Fs)
fprintf("Изменился спектр сигнала: стало больше гармоник (см. графики).\nУ самого сигнала действительно отсеклось все, где амплитуда больше 2, в остальном он остался прежним.\n")
disp('Амплитуда сигнала стала 2, частота осталась прежней.')