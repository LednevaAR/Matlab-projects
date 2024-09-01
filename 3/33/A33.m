clear; clc; close all;
%%
[audio, Fs0] = audioread("task3.wav");
fprintf('Исходная частота дискретизации равна %d\n', Fs0)
figure;
Time = linspace(0, length(audio) ./ Fs0, length(audio));
plot(Time, audio);
xlabel('Time')
ylabel('Amp')
title('Signal')
Spec = fft(audio);
f = (0 : length(audio) - 1) * Fs0 / length(audio);
i = f >= 0 & f <= Fs0/2;
f0 = f(i);
Spec0 = Spec(i);
figure;
plot(f0, abs(Spec0));
xlabel('Frequency')
xlim([0 2000]);
ylabel('Amp')
title('Spec')
%%
Fs = Fs0/2;
%audio1 = zeros(round(size(audio, 1)/2), size(audio, 2));
%for i = 1:length(audio1)
%    audio1(i, :) = audio(2 * i - 1, :);
%end
audio1 = downsample(audio, 2);
audiowrite("1.wav", audio1, Fs);
figure;
Time = linspace(0, length(audio1) ./ Fs, length(audio1));
plot(Time, audio1);
xlabel('Time')
ylabel('Amp')
title('Signal for Fs/2')
Spec = fft(audio1);
f = (0 : length(audio1) - 1) * Fs / length(audio1);
i = f >= 0 & f <= Fs/2;
f1 = f(i);
Spec1 = Spec(i);
figure;
plot(f1, abs(Spec1));
xlabel('Frequency')
ylabel('Amp')
xlim([0 2000]);
title('Spec for Fs/2')
disp('При уменьшении частоты дискретизации качество звука стало хуже, амплитуды гармоник спектра уменьшились.')
%%
Fs = 2 * Fs0;
audio2 = zeros(size(audio, 1) * 2, size(audio, 2));
for i = 1:length(audio2) - 1
    if mod(i, 2) == 0
        audio2(i, :) = (audio(round(i/2), :) + audio(round(i/2) + 1, :))/2;
    else
        audio2(i, :) = audio(round(i/2),:);
    end
end
audiowrite("2.wav", audio2, Fs);
figure;
Time = linspace(0, length(audio2) ./ Fs, length(audio2));
plot(Time, audio2);
xlabel('Time')
ylabel('Amp')
title('Restored signal for 2Fs')
Spec = fft(audio2);
f = (0 : length(audio2) - 1) * Fs / length(audio2);
i = f >= 0 & f <= Fs/2;
f2 = f(i);
Spec2 = Spec(i);
figure;
plot(f2, abs(Spec2));
xlabel('Frequency')
ylabel('Amp')
xlim([0 2000]);
title('Restored spec for 2Fs')
disp('При увеличении частоты дискретизации качество звука стало лучше, амплитуды гармоник спектра увеличились.')
%%
figure;
hold on;
plot(f2, abs(Spec2)/PowerSignal(Spec2));
plot(f1, abs(Spec1)/PowerSignal(Spec1));
plot(f0, abs(Spec0)/PowerSignal(Spec0));
legend('2Fs', 'Fs/2', 'Fs');