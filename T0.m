% Задаем частоту дискретизации и длину сигнала
Fs = 1000; % Герц
t = 0:1/Fs:1-1/Fs;
x = sin(2*pi*50*t) + sin(2*pi*100*t) + sin(2*pi*150*t);

% Уменьшаем частоту дискретизации вдвое
x_downsampled = downsample(x,2);

% Выводим спектр сигнала
X = fft(x);
X_downsampled = fft(x_downsampled);

f = Fs*(0:length(x)-1)/length(x);
f_downsampled = Fs/2*(0:length(x_downsampled)-1)/length(x_downsampled);

figure;
subplot(2,1,1);
plot(f,abs(X));
title('Спектр исходного сигнала');
xlabel('Частота (Гц)');
ylabel('Амплитуда');

subplot(2,1,2);
plot(f_downsampled,abs(X_downsampled));
title('Спектр уменьшенного сигнала');
xlabel('Частота (Гц)');
ylabel('Амплитуда');