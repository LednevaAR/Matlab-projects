clear; clc; close all;
%%
[audio, Fs] = audioread("voice.wav");
Time = linspace(0, length(audio) ./ Fs, length(audio));
figure;
plot(Time, audio);
xlabel('Time')
ylabel('Amp')
title('Signal')
%%
Spec = fft(audio);
Log_Spec = log10(abs(Spec));
f = (0 : length(audio) - 1) * Fs / length(audio);
i = f >= 0 & f <= Fs/2;
f1 = f(i);
Spec1 = Log_Spec(i);
figure;
plot(f1, Spec1);
xlabel('Frequency')
ylabel('Amp')
title('OldLogSpec')
%% этот блок нужен для ликвидации вторичной гармоники, которая была большой и не убиралась весовой функцией
[M, I] = max(Spec);
[M, I] = max(Spec(1 : I - 100));
H = randn(1, 200);
Spec(I - 100 : I + 99) = H / max(H) * (Spec(I - 100) + Spec(I + 99))/2;
%%
audio = ifft(Spec);
audio = audio .* hamming(length(audio));
audio = smoothdata(audio, 'rloess');
%%
Spec = fft(audio);
Log_Spec = log10(abs(Spec));
f = (0 : length(audio) - 1) * Fs / length(audio);
i = f >= 0 & f <= Fs/2;
f1 = f(i);
Spec1 = Log_Spec(i);
figure;
plot(f1, Spec1);
xlabel('Frequency')
ylabel('Amp')
title('NewLogSpec')
[M, I] = max(Spec1);
Ton = (I - 1) * Fs / length(audio);
fprintf('Основной тон моего голосового тракта равен %.0f Гц', Ton);
%%
fig = figure;
f = (0 : length(audio) - 1) * Fs / length(audio);
i = f >= 0 & f <= 2000;
f1 = f(i);
plot(f1, Spec1(i));
xlabel('Frequency')
ylabel('Amp')
xlim([0 2000]);
title('LogSpec')
saveas(fig, "2.jpg")