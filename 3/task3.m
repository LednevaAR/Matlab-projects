clear all; clc; close all;

%% 1 задание
X = ['Задача 1: анализ зашумленного звукового файла'];
disp(X)

[y,Fs] = audioread('file3.wav');

yspec = fft(y);

n = size(y);
yspec0 = fftshift(yspec);         % shift y values
f0 = (-n(1)/2:n(1)/2-1)*(Fs/n(1)); % 0-centered frequency range
power0 = abs(yspec0); 

plot(f0,power0)
xlabel('Frequency')
ylabel('|Y|')
title ("спектр сигнала, задание 1");

% найдем и отфильтруем 3 поврежденные гармоники.

helpMass = yspec0;
yfiltered = zeros(size(yspec0));
for i = 1:3
    max1 = max(abs(helpMass));
    IndexMax1 = find (abs(helpMass) == max1);
    
    X = [' Частота ', num2str(i), ' гармоники: ',num2str(Fs/n(1)*f0(IndexMax1(2)))];
    disp(X)
    yfiltered(IndexMax1) =  helpMass(IndexMax1);
    helpMass(IndexMax1) = 0; 
end

figure
plot(f0,abs(yfiltered))
xlabel('Frequency')
ylabel('|Y_filtered|')
title ("спектр отфильтрованного сигнала, задание 1");

X = ['Отфильтрованный сигнал записан в переменную yfil'];
disp(X)

yfil = ifft(yfiltered);

% sound(yfiltered,Fs);

%% 2 задание
clear all;
X = ['Задача 2: анализ клиппинг эффекта'];
disp(X)

Fs = 2500; % частота дискретизации
T = 3; % продолжительность
N = Fs*T; % количество отсчетов 
F = 1000; % частота в Гц
t = linspace(0,T,N);
amplitude = 3;

y = amplitude*sin(2*pi()*F.*t);

yspec = fft(y);

n = size(y);
yspec0 = fftshift(yspec);         % shift y values
f0 = (-n(2)/2:n(2)/2-1)*(Fs/n(2)); % 0-centered frequency range
power0 = abs(yspec0); 

figure
plot(f0,power0)
xlabel('Frequency')
ylabel('|Y|')
title ("спектр созданного сигнала, задание 2");

X = ['Создадим звуковой файл Sin.wav и запишем туда созданный сигнал'];
disp(X)

audiowrite("Sin.wav",y./amplitude,Fs)

% Сымитируем клипинг эффект и запишем результат в матрицу y

%   sound(y,Fs);
Umax = 1; 

y(y>Umax) = Umax;
y(y<-Umax) = -Umax;

yspecClipping = fft(y);

n = size(y);
yspec0 = fftshift(yspecClipping);         % shift y values
f0 = (-n(2)/2:n(2)/2-1)*(Fs/n(2)); % 0-centered frequency range
power0 = abs(yspec0); 

figure
plot(f0,power0)
xlabel('Frequency')
ylabel('|Y|')
title ("спектр сигнала с клипинг эффектом, задание 2");


X = ['В частотном спектре клиппирование приводит к увеличению количества гармоник в области низких частот'];
disp(X)
X = ['Сам сигнал имеет большую площадь пиков, чем максимальный не клиппированный'];
disp(X)

%   sound(y,Fs);

%% Задача 3: анализ влияния частоты дискретизации.
clear all;

X = ['Задача 3: анализ влияния частоты дискретизации'];
disp(X)

[yTask3,Fs] = audioread('task3.wav');

X = [' Частота дискретизации ', num2str(Fs)];
disp(X)

yspecClipping = fft(yTask3);

% уменьшим частоту дискретизацию вдвое 

yTask3Discr = zeros(fix(length(yTask3')/2),2);
yTask3Discr(1:end, 1:end) =  yTask3(2:2:end, 1:end); 

X = ['Уменьшения частоты дискретизации приводит к заметному ухудшению качества звука, спекрт сигнала зашумлен'];
disp(X)

% sound(yTask3Discr,Fs3/2);

% увеличим дискретизацию вдвое 

yTask3Max = zeros((length(yTask3')*2),2);

X = ['Увеличение частоты дискретизации улучшает качество звука, спекрт сигнала содержит меньше гармоник'];
disp(X)

for i = 1:length(yTask3')-1
    yTask3Max(2*i-1,1:end) = yTask3(i,1:end);
    yTask3Max(2*i,1:end) = (yTask3(i,1:end) + yTask3(i+1,1:end))./2;
end

% Построение графиков
figure

yTask3Power = PowerSignal(yspecClipping);
n = size(yTask3);
yspec0 = fftshift(yspecClipping(1:end,1))./yTask3Power(1);         % shift y values
f0 = (-n(1)/2:n(1)/2-1)*(Fs/n(1)); % 0-centered frequency range
power0 = abs(yspec0); 
h = fix(length(yspec0)/2);

plot(f0(1,h:end),power0(h:end,1))


hold on


Fs = Fs/2;
yTask3DiscrSpec = fft(yTask3Discr);
yTask3DiscrPower = PowerSignal(yTask3DiscrSpec);

n = size(yTask3Discr);
yspec0 = fftshift(yTask3DiscrSpec(1:end,1))./yTask3DiscrPower(1);         % shift y values
f0 = (-n(1)/2:n(1)/2-1)*(Fs/n(1)); % 0-centered frequency range
power0 = abs(yspec0);
h = fix(length(yspec0)/2);

plot(f0(1,h:end),power0(h:end,1))


hold on

Fs = Fs*4;
yTask3MaxSpec = fft(yTask3Max);
yTask3MaxPower = PowerSignal(yTask3MaxSpec);

n = size(yTask3Max);
yspec0 = fftshift(yTask3MaxSpec(1:end,1))./yTask3MaxPower(1);         % shift y values
f0 = (-n(1)/2:n(1)/2-1)*(Fs/n(1)); % 0-centered frequency range
power0 = abs(yspec0); 
h = fix(length(yspec0)/2);

plot(f0(1,h:end),power0(h:end,1))

xlabel('Frequency')
ylabel('|Y|/Power')
title ("спектр сигналов, задание 3");
legend ("исходный сигнал","частота дискретизации уменьшенна вдвое","частота дискретизации увеличена вдвое");
hold off
% sound(yTask3Max,Fs3*2);

function power = PowerSignal(Signal)
    power = mean(Signal.*conj(Signal));
end