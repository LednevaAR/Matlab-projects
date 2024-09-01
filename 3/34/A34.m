clear; clc; close all;
%%
[audio, Fs] = audioread("budgie-chirping.wav");
%figure
%spectrogram(audio(:, 1), 'yaxis')
%figure
%spectrogram(audio(:, 2), 'yaxis')
figure
Time = linspace(0, length(audio) ./ Fs, length(audio));
plot(Time, audio);
xlabel('Time')
ylabel('Amp')
title('Signal')
n = length(audio);
f = (-n/2:n/2-1)*(Fs/n);
%Spec = fft(audio);
%figure
%plot(f, abs(Spec));
%xlabel('Frequency')
%ylabel('Amp')
%title('Spec')
%%
result = audio;
i = find(abs(audio(:, 1)) > 0.6);
j = find(abs(audio(:, 2)) > 0.6);
noise = zeros(size(audio));
x = linspace(0, Fs, n);
for k = 1:size(i, 1)
    noise(i(k) - 100 : i(k) + 100, 1) = randn(201, 1);
    result(:, 1) = result(:, 1) + noise(:, 1) * audio(i(k), 1)/max(abs(noise(:, 1)));
    noise = zeros(size(audio));
end
for k = 1:size(j, 1)
    noise(j(k) - 100 : j(k) + 100, 1) = randn(201, 1);
    result(:, 2) = result(:, 2) + noise(:, 2) * audio(j(k), 2)/max(abs(noise(:, 2)));
    noise = zeros(size(audio));
end
result(:, 1) = result(:, 1)/max(abs(result(:, 1)));
result(:, 2) = result(:, 2)/max(abs(result(:, 2)));
audiowrite('result.wav', result, Fs);