close all; clear; clc;
%% Init parametrs of model
Length_Bit_vector = 12e5;

%Constellation = "BPSK"; % QPSK, 8PSK, 16-QAM
%Constellation = "QPSK";
Constellation = "8PSK";
%Constellation = "16-QAM";

SNR = 30; % dB

%% Bit generator

Bit_Tx = randi([0,1], 1, Length_Bit_vector);
%20*log10(mean(abs(Tx_OFDM_Signal))/mean(abs(Rx_OFDM_Signal - Tx_OFDM_Signal)))
IQ_TX = mapping(Bit_Tx, "16-QAM");
IQ_RX = NoiseGenerator(SNR, IQ_TX);
MER_estimation = MER_my_func(IQ_RX, "16-QAM");
%tic;
%transmit(Bit_Tx);
%Time = toc;
%T = Time/Length_Bit_vector;
%% Mapping

IQ_TX = mapping(Bit_Tx, Constellation);
figure();
scatter(real(IQ_TX), imag(IQ_TX), "filled")
title("Scatter plot")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, Bit_depth_Dict] = constellation_func(Constellation);
Dict = name(Constellation);
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
%% Channel
% Write your own function Eb_N0_convert(), which convert SNR to Eb/N0
Eb_N0 = Eb_N0_convert(SNR, Constellation);%
%S = convertSNR(Eb_N0, "ebno", BitsPerSymbol = 3);
% Use your own function of generating of AWGN from previous tasks
IQ_RX = NoiseGenerator(SNR, IQ_TX);
figure();
scatter(real(IQ_RX), imag(IQ_RX), "filled")
title("Scatter plot")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
%% Demapping
Bit = demapping(IQ_RX, Constellation);
Bit_Rx = Bit(1:length(Bit_Tx));
%% Error check
% Write your own function Error_check() for calculation of BER
BER = Error_check(Bit_Tx, Bit_Rx);
fprintf("Вероятность ошибки равна BER = %f", BER)
%% Additional task. Modulation error ration
SNR0 = -50:1:50;
C = ["BPSK", "QPSK", "8PSK", "16-QAM"];
MER0 = zeros(4, length(SNR0));
for j = 1:4 
    for i = 1:length(SNR0)
        Bit_Tx = randi([0,1], 1, Length_Bit_vector);
        IQ_TX = mapping(Bit_Tx, C(j));
        IQ_RX = NoiseGenerator(SNR0(i), IQ_TX);
        MER_estimation = MER_my_func(IQ_RX, C(j));
        MER0(j, i) = MER_estimation;
    end
end
f = figure;
plot(SNR0, abs(MER0(1, 1:end) - SNR0), 'DisplayName','BPSK')
hold on
plot(SNR0, abs(MER0(2, 1:end) - SNR0), 'DisplayName','QPSK')
plot(SNR0, abs(MER0(3, 1:end) - SNR0), 'DisplayName','8PSK')
plot(SNR0, abs(MER0(4, 1:end) - SNR0), 'DisplayName','16QAM')
hold off
grid on
title("MER")
legend
ylabel('MER (dB)')
xlabel('SNR (dB)')
saveas(f, "MER.fig");
% Compare the SNR and MER_estimation from -50dB to +50dB for BPSK, QPSK,
% 8PSK and 16QAM constellation.
% Plot the function of error between SNR and MER for each constellation 
% Discribe the results. Make an conclusion about MER.
% You can use the cycle for collecting of data
% Save figure
%% Experimental BER(SNR) and BER(Eb/N0)
% Collect enough data to plot BER(SNR) and BER(Eb/N0) for each
% constellation.
% Compare the constalation. Describe the results
% You can use the cycle for collecting of data
% Save figure
SNR0 = 0:0.1:22;
BER0 = zeros(4, length(length(SNR0)));
EbN00 = zeros(4, length(length(SNR0)));
C = ["BPSK", "QPSK", "8PSK", "16-QAM"];
for j = 1:4
    for i = 1:length(SNR0)
        Bit_Tx = randi([0,1], 1, Length_Bit_vector);
        EbN00(j, i) = Eb_N0_convert(SNR0(i), C(j));    
        IQ_TX = mapping(Bit_Tx, C(j));
        IQ_RX = NoiseGenerator(SNR0(i), IQ_TX);
        Bit = demapping(IQ_RX, C(j));
        Bit_Rx = Bit(1:length(Bit_Tx));
        BER = Error_check(Bit_Tx, Bit_Rx);
        BER0(j, i) = BER;
    end
end
f = figure;
semilogy(SNR0, BER0(1, 1:end), 'DisplayName','BPSK')
hold on
semilogy(SNR0, BER0(2, 1:end), 'DisplayName','QPSK')
semilogy(SNR0, BER0(3, 1:end), 'DisplayName','8PSK')
semilogy(SNR0, BER0(4, 1:end), 'DisplayName','16QAM')
hold off
grid on
title("Experimental lines of BER(SNR)")
legend
ylabel('BER')
xlabel('SNR (dB)')
ylim([10^(-5) 10^0])
xlim([0 18])
saveas(f, "ExperimentSNR.fig")
f = figure;
semilogy(EbN00(1, 1:end), BER0(1, 1:end), 'DisplayName','BPSK')
hold on
semilogy(EbN00(2, 1:end), BER0(2, 1:end), 'DisplayName','QPSK')
semilogy(EbN00(3, 1:end), BER0(3, 1:end), 'DisplayName','8PSK')
semilogy(EbN00(4, 1:end), BER0(4, 1:end), 'DisplayName','16QAM')
hold off
grid on
title("Experimental lines of BER(Eb/N0)")
legend
ylabel('BER')
xlabel('E_b/N_0 (dB)')
ylim([10^(-8) 10^0])
xlim([0 18])
saveas(f, "Experiment.fig")
%% Theoretical lines of BER(Eb/N0)
% Read about function erfc(x) or similar
% Configure the function and get the theoretical lines of BER(Eb/N0)
% Compare the experimental BER(Eb/N0) and theoretical for BPSK, QPSK, 8PSK
% and 16QAM constellation
% Save figure
EbN0_dB = 0:0.1:22;
EbN0 = 10 .^ (EbN0_dB/10);
BER2 = (1/2) * erfc(sqrt(EbN0));
BER4 = (1/2) * erfc(sqrt(EbN0));
BER8 = (1/3) * erfc(sqrt(EbN0 * 3) * sin(pi/8));
BER16 = (1/2) * erfc(sqrt(EbN0 * 6 / 15));
f = figure;
semilogy(EbN0_dB, BER2, 'DisplayName','BPSK/QPSK')
hold on
semilogy(EbN0_dB, BER4, 'DisplayName','BPSK/QPSK')
semilogy(EbN0_dB, BER8, 'DisplayName','8PSK')
semilogy(EbN0_dB, BER16, 'DisplayName','16QAM')
hold off
grid on
title("Theoretical lines of BER(Eb/N0)")
legend
ylabel('BER')
xlabel('E_b/N_0 (dB)')
ylim([10^(-8) 10^0])
xlim([0 18])
saveas(f, "Theory.fig")
%%
f = figure;
semilogy(EbN0_dB, BER2, 'DisplayName','Theory')
hold on
semilogy(EbN00(1, 1:end), BER0(1, 1:end), 'DisplayName','Experiment')
hold off
grid on
title("BPSK")
legend
ylabel('BER')
xlabel('E_b/N_0 (dB)')
ylim([10^(-8) 10^0])
xlim([0 18])
saveas(f, "Theory1.fig")
%%
f = figure;
semilogy(EbN0_dB, BER4, 'DisplayName','Theory')
hold on
semilogy(EbN00(2, 1:end), BER0(2, 1:end), 'DisplayName','Experiment')
hold off
grid on
title("QPSK")
legend
ylabel('BER')
xlabel('E_b/N_0 (dB)')
ylim([10^(-8) 10^0])
xlim([0 18])
saveas(f, "Theory2.fig")
%%
f = figure;
semilogy(EbN0_dB, BER8, 'DisplayName','Theory')
hold on
semilogy(EbN00(3, 1:end), BER0(3, 1:end), 'DisplayName','Experiment')
hold off
grid on
title("8PSK")
legend
ylabel('BER')
xlabel('E_b/N_0 (dB)')
ylim([10^(-8) 10^0])
xlim([0 18])
saveas(f, "Theory3.fig")
%%
f = figure;
semilogy(EbN0_dB, BER16, 'DisplayName','Theory')
hold on
semilogy(EbN00(4, 1:end), BER0(4, 1:end), 'DisplayName','Experiment')
hold off
grid on
title("16-QAM")
legend
ylabel('BER')
xlabel('E_b/N_0 (dB)')
ylim([10^(-8) 10^0])
xlim([0 18])
saveas(f, "Theory4.fig")