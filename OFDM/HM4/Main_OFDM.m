clc; % чистка командного окна
close all; % закрыть дополнительные окна 
clear all; % очистить память
rng(1); % фиксирование начального состояния генератора случайных чисел Матлаба
%% 
% Параметры
% Конфигурация модели
% Выбор созвездия
N_carrier = 400;
N_fft = 1024;
T_guard = N_fft / 8;
Amount_OFDM_Frames = 60;
Amount_OFDM_Symbols_per_Frame = 5;
Time_Delay = randi(N_fft + T_guard + 1) - 1;
%Freq_Shift = (randi(61) - 31);
Freq_Shift = 0.1;
Channel = [0, 1; 4, 0.6; 10, 0.3];
SNR = randi(101) - 1;
% SNR = zeros(1, 31);
Percent_pilot = 20;
amp_pilots = 4/3;
%%
PAPR = zeros(1, 3);
CCDF = zeros(3, (Amount_OFDM_Symbols_per_Frame * Amount_OFDM_Frames - 1) * (N_fft + T_guard) + 1);
PAPR_sorted = zeros(3, (Amount_OFDM_Symbols_per_Frame * Amount_OFDM_Frames - 1) * (N_fft + T_guard) + 1);
n = 1;
f = figure;
for phase_mod = [1 2 3]
%for Percent_pilot = [0 5 10 20 25]
if Percent_pilot == 0
    pilot_index = [];
else
    pilot_index = [1 : floor(100 / Percent_pilot) : N_carrier - floor(100 / Percent_pilot) / 2, N_carrier];
end
inform_index = setdiff(1 : 1 : N_carrier, pilot_index);
% Количество пилотов в полосе
amount_ration_pilots = length(pilot_index);
%%
if_CCDF = true;
%if if_CCDF == true
    %f = figure;
    var_randomizer = true;%[false, true];
% else
%     var_randomizer = true;
% end
%%
BER = zeros(size(SNR));
% for j = 1 : 31
% SNR(j) = j - 1;
for if_randomizer = var_randomizer
Register = [1 0 0 1 0 1 0 1 0 0 0 0 0 0 0];
constellation = "16-QAM"; % "BPSK"; % "QPSK"
File = 'HM1.jpg'; % Адрес файла
[Dictionary, D, Bit_depth_Dict] = constellation_func(constellation); % to-do lab 1
QAM_cells = length(Dictionary); % количество точек созвездия
%%
%Передатчик
len = Amount_OFDM_Frames * Amount_OFDM_Symbols_per_Frame * length(inform_index) * Bit_depth_Dict;
%Чтение из файлы последовательности бит длиной len
Input_Bit_Buffer = file_reader(File, len);
%disp(sum(Input_Bit_Buffer > 255))
%Рандомизация
if if_randomizer == true
    Input_Bit_Buffer_if_randomizer = randomizer(Input_Bit_Buffer, Register, Amount_OFDM_Frames);
else
    Input_Bit_Buffer_if_randomizer = Input_Bit_Buffer;
end
%Маппинг
Tx_IQ_points = mapping(Input_Bit_Buffer_if_randomizer, constellation);
%Создание OFDM символов
Tx_OFDM_symbols = OFDM_Mod(Tx_IQ_points, N_fft, T_guard, inform_index, pilot_index, amp_pilots, phase_mod);
%Сшивка OFDM символов
Tx_OFDM_Signal = signal_generator(Tx_OFDM_symbols);
%%
%Расчет PAPR и CCDF, построение графика CCDF(PAPR)
if if_CCDF == true
    PAPR(n) = 10 * log10(max(abs(Tx_OFDM_Signal) .^ 2) / mean(abs(Tx_OFDM_Signal) .^ 2));
    PAPR_slid_wind = compute_PAPR_slid_wind(Tx_OFDM_Signal, N_fft + T_guard);
    [CCDF(n, :), PAPR_sorted(n, :)] = compute_CCDF(PAPR_slid_wind);
    semilogy(PAPR_sorted(n, :), CCDF(n, :));
    n = n + 1;
    hold on
end
end
end
%%
% f = figure;
% semilogy(PAPR_sorted(1, :), CCDF(1, :));
% hold on
% semilogy(PAPR_sorted(2, :), CCDF(2, :));
% hold on
% semilogy(PAPR_sorted(3, :), CCDF(3, :));
% hold on
% semilogy(PAPR_sorted(5, :), CCDF(5, :));
if if_CCDF == true
title("CCDF(PAPR)")
xlabel("PAPR, dB")
ylabel("CCDF")
saveas(f, "CCDF_PAPR.fig")
end
%%
%канал
Rx_OFDM_Signal = Tx_OFDM_Signal;
% Гауссовский шум
%Rx_OFDM_Signal = NoiseGenerator(SNR, Rx_OFDM_Signal);
% Частотная рассинхронизация
%Rx_OFDM_Signal = Rx_OFDM_Signal .* exp(1j * 2 * (1 : size(Rx_OFDM_Signal, 2)) * pi * Freq_Shift  / N_fft);
% Многолучевое распространение
%Rx_OFDM_Signal = multipath(Rx_OFDM_Signal, Channel);
% Временная рассинхронизация
%Rx_OFDM_Signal = [Rx_OFDM_Signal(1, Time_Delay + 1 : end), zeros(1, Time_Delay)];
%%
%приемник
%Rx_OFDM_Signal = [Rx_OFDM_Signal(1 + T_guard : end), zeros(1, T_guard)]; %T_guard / 2, T_guard
Rx_OFDM_data = OFDM_Signal_Demod(Rx_OFDM_Signal, T_guard, N_fft);
Rx_IQ = zeros(size(Rx_OFDM_data, 1), N_fft);
for i = 1 : size(Rx_OFDM_data, 1)
    Rx_IQ(i, 1 : N_fft) = fft(Rx_OFDM_data(i, 1 : end), N_fft);
end
Rx_IQ_points = conj(reshape(Rx_IQ(:, 1 : N_carrier)', 1, numel(Rx_IQ(:, 1 : N_carrier))));
%Rx_IQ_points = conj(reshape(Rx_IQ(:, inform_index)', 1, numel(Rx_IQ(:, inform_index))));
%Вычисление MER
%MER = MER_my_func(Rx_IQ_points, constellation);
%Демаппинг
Output_Bit_Buffer_if_randomizer = demapping(Rx_IQ_points, constellation);
%Дерандомизация
if if_randomizer == true
    Output_Bit_Buffer = randomizer(Output_Bit_Buffer_if_randomizer, Register, Amount_OFDM_Frames);
else
    Output_Bit_Buffer = Output_Bit_Buffer_if_randomizer;
end
%Вычисление, какая часть бит была получена верно
%Probability = Error_check(Input_Bit_Buffer, Output_Bit_Buffer);
%Вычисление BER
% BER(j) = 1 - Error_check(Input_Bit_Buffer, Output_Bit_Buffer);
% end
%%
%Построение графика BER(SNR)
% f = figure();
% semilogy(SNR, BER)
% title("BER(SNR)")
% xlabel("SNR, dB")
% ylabel("BER")
% grid on
% saveas(f, "BER_SNR.fig")
%%
%Построение графика АЧХ
f = figure();
plot(1 : size(Rx_IQ, 2), abs(Rx_IQ(2, 1 : end)))
xlim([0 1024])
title("Plot")
xlabel("Frequency")
ylabel("Amplitude")
grid on
saveas(f, "АЧХ.fig")
%Построение сигнального созвездия
f = figure();
s = scatter(real(Rx_IQ_points), imag(Rx_IQ_points), "filled");
s.SizeData = 5;
title("Plot")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, Bit_depth_Dict] = constellation_func(constellation);
Dict = name(constellation);
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
grid on
saveas(f, "Const.fig")
%% 
%Восстановление картинки
% Output_Buffer = num2str(reshape(Output_Bit_Buffer, numel(Output_Bit_Buffer) / 8, 8));
% for i = 1 : size(Output_Buffer, 1)
%     Output(i, 1 : 8) = regexprep(Output_Buffer(i, 1 : end), ' ', '');
% end
% Output = bin2dec(Output);
% Output = Output';
% Out = uint8(reshape(Output, 100, 200, 3));
% imwrite(Out, 'image.jpg','jpg');