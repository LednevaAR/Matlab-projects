% ����������.
% =========================================================================
%> ���������� �������� �����
% =========================================================================
    %> �������� workspace
    clear all;
    %> �������� ��������
    close all;
    %> �������� Command Window
    clc;
    %%
% =========================================================================
%> ������� sinc. ������ �� ������.
% =========================================================================
    %> ������� ������
    x = [-10:0.1:10];
    %> ������� sinc
    y = sinc(x);
    % =====================================================================
    %> ������ ���������� ��������������
    % =====================================================================
    f = figure;
    plot(x,y)
    title('sinc(x)')
    saveas(f, "Part_0_impulse.fig")
    % =====================================================================
    %> ������ �������
    % =====================================================================
    %> ����� + ������ � ��
    spectum = 10*log10(abs(fft(y)));
    %> ������ �� ������
    spectum = [spectum(102:201), spectum(1:101)];
    %> ������ � ������ T = 2 (��. �������� sinc)
    figure 
    plot(x/2, spectum)
    title('spectum sinc(x)')
    saveas(f, "Part_0_spectrum.fig")
%%
% =========================================================================
%> ������ 1: �������� �������, ������� ���������� ����������� (���������� 
%> ��������������)��� ������� ������ �� ������������ ��������.
%> ��������� ���������� � ��������� �������������� �������.
% =========================================================================
    %> ����� ������� � �������� (����� ������� ��������� sinc, ����� � ���� ������)
    span = 20;
    %> ����� ������� �� ������
    nsamp  = 4;
    %> ���������� ����������� (alfa)
    rolloff = 0.2;
    % =====================================================================
    %> @todo ��������� �������
    sqimpuls = sqRCcoeff (span, nsamp, rolloff);
    %> @todo ��������� ���������� � ��������� �������������� �������
    % .........
    x = -span/2:(1/nsamp):span/2;
    f = figure;
    plot(x, sqimpuls)
    title('Task 1, impuls')
    saveas(f, "Part_1_impuls.fig")
    spectrum = 10*log10(abs(fft(sqimpuls)));
    spectrum = [spectrum((span * nsamp / 2 + 2):(span * nsamp + 1)), spectrum(1:(span * nsamp / 2 + 1))];
    f = figure;
    plot(x, spectrum)
    title('Task 1, spectrum, dB')
    saveas(f, "Part_1_spectrum.fig")
% =========================================================================
%> �������� 1.
%> ��������� �� ����������� ��������
% =========================================================================
txfilter1 = comm.RaisedCosineTransmitFilter('RolloffFactor', rolloff, ...
                                           'FilterSpanInSymbols',span,...
                                           'OutputSamplesPerSymbol', nsamp);
check1 = coeffs(txfilter1);
if sum(abs(check1.Numerator-sqimpuls))< 0.001 % �������� ���������� ���� 
                                              % ���������� ������������� 
                                              % � �������� ���������
    ans = '�������� ������ 1 �������� �������'
else 
    err = '������ � ������ 1. ��������� �����������'
    ans = sum(abs(check1.Numerator-sqimpuls))
end
%%
% =========================================================================
%> ������ 2: �������� �������, ������� ���������� ����������� (���������� 
%> ��������������) ��� ������� ������������ ��������.
%> ��������� ���������� � ��������� �������������� �������.
%> ��������� ���������� �������������� ��� �����, ��� ����� � �������������� sinc 
%> �� ����� �������
% =========================================================================
    %> ����� ������� � �������� (����� ������� ��������� sinc, ����� � ���� ������)
    span = 20;
    %> ����� ������� �� ������
    nsamp  = 4;
    %> ���������� ����������� (alfa)
    rolloff = 0.2;
    % =====================================================================
    %> @todo ��������� �������
    impuls = RCcoeff (span, nsamp, rolloff);
    %> @todo ��������� ���������� � ��������� �������������� �������
    % .........
    x = -span/2:(1 / nsamp):span/2;
    f = figure;
    plot(x, impuls)
    title('Task 2, impuls')
    saveas(f, "Part_2_impuls.fig")
    spectrum = 10*log10(abs(fft(impuls)));
    spectrum = [spectrum((span * nsamp / 2 + 2):(span * nsamp + 1)), spectrum(1:(span * nsamp / 2 + 1))];
    f = figure;
    plot(x, spectrum)
    title('Task 2, spectrum, dB')
    saveas(f, "Part_2_spectrum.fig")
    
% =========================================================================
%> �������� 2.
%> ��������� �� ����������� ��������
% =========================================================================
txfilter2 = comm.RaisedCosineTransmitFilter('RolloffFactor', rolloff, ...
                                            'FilterSpanInSymbols',span,...
                                            'OutputSamplesPerSymbol', nsamp,...
                                            'Shape', 'Normal');
check2 = coeffs(txfilter2);
if sum(abs(check2.Numerator-impuls))< 0.1 % �������� ���������� ���� 
                                          % ���������� ������������� 
                                          % � �������� ���������
    ans = '�������� ������ 2 �������� �������'
else 
    err = '������ � ������ 2. ��������� �����������'
    ans = sum(abs(check2.Numerator-impuls))
end
%%
% =========================================================================
%> ������� 3. 
%> �������� ������� ����������, ������� �������� � ���� �������: �
%> ����������� ���������� ������� �� ������ � ��� (��������� ����������)
%> @warning ������������ ������� mapping �� ������� �����
% =========================================================================
    UpSempFlag = true(1);
    bits = randi([0 1], 1, 1000); % ��������� ���
    sign = mapping (bits, "QPSK");       %QPSK 500 �������� 
    filtsign = Filtration(sign, sqimpuls, nsamp, UpSempFlag);
    % =====================================================================
    %> �������� 3.1
    %> ��������� ������������ ������ �-�� � ������������������ �� ������������ ��������.
    % =====================================================================
    check3 = txfilter1(sign.').';
    if sum(abs(check3-filtsign))< 0.1 % �������� ���������� ���� 
                                      % ���������� ������������� 
                                      % � �������� ���������
        ans = '�������� ������ 3.1 �������� �������'
    else 
        err = '������ � ������ 3.1. ��������� ������'
        ans = sum(abs(check3-filtsign))
    end
    
    % =====================================================================
    %> �������� 3.2
    %> ��������� ������������ ������ �-�� ��� ����������������� �� ������������ ��������.
    % =====================================================================
    UpSempFlag = false(1);
    filtsign2 = Filtration(filtsign, sqimpuls, nsamp, UpSempFlag);
    rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor', rolloff, ...
                                              'FilterSpanInSymbols',span,...
                                              'InputSamplesPerSymbol', nsamp,...
                                              'DecimationFactor', 1);
    check4 = rxfilter(filtsign.').';
    if sum(abs(check4-filtsign2))< 0.1 % �������� ���������� ���� 
                                       % ���������� ������������� 
                                       % � �������� ���������
        ans = '�������� ������ 3.2 �������� �������'
    else 
        err = '������ � ������ 3.2. ��������� ������'
        ans = sum(abs(check4-filtsign2))
    end
%% 3
nsamp  = 4;
impuls = RCcoeff (span, nsamp, rolloff);
bits = randi([0 1], 1, 1000);
sign = mapping(bits, "QPSK");
UpSempFlag = true(1);
filtsign31 = Filtration(sign, impuls, nsamp, UpSempFlag);
f = figure();
scatter(real(filtsign31), imag(filtsign31), "filled")
title("Plot before")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, Bit_depth_Dict] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
saveas(f, "Const_3_before.fig")
f = figure();
plot(abs(filtsign31))
title("Plot before")
xlabel("Time")
ylabel("Ampl")
ylim([0 5])
grid on
saveas(f, "Ampl_3_before.fig")

UpSempFlag = false(1);
filtsign32 = Filtration(filtsign31, impuls, nsamp, UpSempFlag);
f = figure();
scatter(real(filtsign32), imag(filtsign32), "filled")
title("Plot after")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, Bit_depth_Dict] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
saveas(f, "Const_3_after.fig")
f = figure();
plot(abs(filtsign32))
title("Plot after")
xlabel("Time")
ylabel("Ampl")
ylim([0 5])
grid on
saveas(f, "Ampl_3_after.fig")
%% 4
bits = randi([0 1], 1, 1000);
sign = mapping(bits, "QPSK");
UpSempFlag = true(1);
filtsign41 = Filtration(sign, impuls, nsamp, UpSempFlag);
f = figure();
scatter(real(filtsign41), imag(filtsign41), "filled")
title("Plot before")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, Bit_depth_Dict] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
saveas(f, "Const_4_before.fig")
f = figure();
plot(abs(filtsign41))
title("Plot before")
xlabel("Time")
ylabel("Ampl")
ylim([0 5])
grid on
saveas(f, "Ampl_4_before.fig")

%filtsign41 = NoiseGenerator(10, filtsign41);

UpSempFlag = false(1);
filtsign42 = Filtration(filtsign41, impuls, nsamp, UpSempFlag);
filtsign42 = filtsign42(1 : nsamp : end);
%filtsign42 = resample(filtsign42, 1, 4);
f = figure();
scatter(real(filtsign42), imag(filtsign42), "filled")
title("Plot after")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, Bit_depth_Dict] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
saveas(f, "Const_4_after.fig")
f = figure();
plot(abs(filtsign42))
title("Plot after")
xlabel("Time")
ylabel("Ampl")
ylim([0 5])
grid on
saveas(f, "Ampl_4_after.fig")
%% 5
bits = randi([0 1], 1, 1000);
sign = mapping(bits, "QPSK");
UpSempFlag = true(1);
filtsign51 = Filtration(sign, sqimpuls, nsamp, UpSempFlag);
SNR = 30;
filtsign51 = NoiseGenerator(SNR, filtsign51);
freq_offset_percent = linspace(-150, 150, 301);
MER = zeros(1, length(freq_offset_percent));
time = (0 : length(filtsign51) - 1) / nsamp;
for i = 1 : length(freq_offset_percent)
    filtsign5 = filtsign51 .* exp(1j * 2 * pi * time * freq_offset_percent(i) / 100);
    UpSempFlag = false(1);
    filtsign52 = Filtration(filtsign5, sqimpuls, nsamp, UpSempFlag);
    filtsign52 = filtsign52 .* exp(-1j * 2 * pi * time * freq_offset_percent(i) / 100);
    filtsign52 = filtsign52(1 : nsamp : end);
    %filtsign52 = resample(filtsign52, 1, nsamp);
    MER(i) = MER_my_func(filtsign52, "QPSK");
end
f = figure();
plot(freq_offset_percent, MER - SNR)
xlabel("Freq offset percent")
ylabel("MER - SNR")
grid on
saveas(f, "MER_5.fig")
%%
bits = randi([0 1], 1, 1000);
sign = mapping(bits, "QPSK");
UpSempFlag = true(1);
filtsign51 = Filtration(sign, sqimpuls, nsamp, UpSempFlag);
SNR = 30;
filtsign51 = NoiseGenerator(SNR, filtsign51);
freq_offset_percent = 20; %linspace(-150, 150, 301);
time = (0 : length(filtsign51) - 1) / nsamp;
filtsign5 = filtsign51 .* exp(1j * 2 * pi * time * freq_offset_percent / 100);
UpSempFlag = false(1);
filtsign52 = Filtration(filtsign5, sqimpuls, nsamp, UpSempFlag);
filtsign52 = filtsign52 .* exp(-1j * 2 * pi * time * freq_offset_percent / 100);
filtsign52 = filtsign52(1 : nsamp : end);
%filtsign52 = resample(filtsign52, 1, nsamp);
f = figure();
scatter(real(filtsign52), imag(filtsign52), "filled")
title("Plot after")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, Bit_depth_Dict] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
%%
bits = randi([0 1], 1, 1000);
sign = mapping(bits, "QPSK");
UpSempFlag = true(1);
filtsign41 = Filtration(sign, sqimpuls, nsamp, UpSempFlag);
%filtsign41 = NoiseGenerator(30, filtsign41);

UpSempFlag = false(1);
filtsign42 = Filtration(filtsign41, sqimpuls, nsamp, UpSempFlag);
filtsign42 = filtsign42(2 : nsamp : end);
%filtsign42 = resample(filtsign42, 1, 4);
f = figure();
plot(abs(filtsign42))
title("Plot after")
xlabel("Time")
ylabel("Ampl")
grid on
f = figure();
scatter(real(filtsign42), imag(filtsign42), "filled")
title("Plot after")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, Bit_depth_Dict] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on