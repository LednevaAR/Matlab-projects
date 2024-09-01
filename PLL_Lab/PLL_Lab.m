%%
close; clear; clc;

%% config
Freq_Offset = 0.1; % normalised frequency
SNR = 30; % dB
%% Transmitter
% config
Amount_of_Frame = 1000;
Length_Data_IQ = 1440;

% Start of frame
SOF = [1 0 0 1 1 1 0 1 0 1 0 1 0 1 1 0 0 1 0 0]; 
IQ_SOF = mapping(SOF, "BPSK"); % Use this sequence on the Rx as a Pilot-Signal

% bit generator


% QAM | mapper
% Generate vector of bits with length = Amount_of_Frame*Length_Data_IQ
Tx_Bits = randi([0 1], Amount_of_Frame * Length_Data_IQ * 2, 1)'; 
% Use your function Mapping.m from previous courses
TX_IQ_Data = mapping(Tx_Bits, "QPSK");



% Frame structure 
% |20 IQ BPSK Start-of-Frame| 1440 IQ QPSK ... 
IQ_TX_Frame = FrameStruct(TX_IQ_Data, IQ_SOF, Amount_of_Frame);



%% Channel
% Add white Gaussian noise to signal
%Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
% Add frequency offset
%Channel_IQ1 = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2))* pi * Freq_Offset);
%%
% Configurate the Loop Filter
% =========================================================================
% Loop filter preset
% -------------------------------------------------------------------------
    Xi = 1;           % detector gainDampingFactor
    BnTs = 0.5;       % Normalized loop bandwidth (Ts = 1 for this problem)
    Kd = 2*pi;          % Phase (not change)
    K0 = 1;             % not change
% =========================================================================
%> Loop filter coefficient calculation
% -------------------------------------------------------------------------
    %todo прописать расчет коэффициентов фильтра
    %> Proportional coefficient
    Kp = 2 * Xi * (BnTs / (Xi + 1 / (4 * Xi))) / (Kd * K0);
    %> Integrator coefficient
    Ki = (BnTs / (Xi + 1 / (4 * Xi))) ^ 2 / (Kd * K0);

Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * Freq_Offset);
[RX_IQ_DM, DM_estimate, DM_filtered] = DM(Channel_IQ, Kp, Ki, IQ_SOF);
f = figure();
scatter(real(Channel_IQ(1000, 1 : end)), imag(Channel_IQ(1000, 1 : end)), "filled")
title("Plot before")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, ~] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
saveas(f, "Const_DM_before.fig")
f = figure();
scatter(real(RX_IQ_DM(1000, 1 : end)), imag(RX_IQ_DM(1000, 1 : end)), "filled")
title("Plot after")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, ~] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
saveas(f, "Const_DM_after.fig")


%% Receiver with frequency estimator based on Delay and Multiply with D=2
% in feedback scheme

% Configurate the Loop Filter
% =========================================================================
% Loop filter preset
% -------------------------------------------------------------------------
    Xi = 1;           % detector gainDampingFactor
    BnTs = 0.5;       % Normalized loop bandwidth (Ts = 1 for this problem)
    Kd = 2*pi;          % Phase (not change)
    K0 = 1;             % not change
% =========================================================================
%> Loop filter coefficient calculation
% -------------------------------------------------------------------------
    %todo прописать расчет коэффициентов фильтра
    %> Proportional coefficient
    Kp = 2 * Xi * (BnTs / (Xi + 1 / (4 * Xi))) / (Kd * K0);
    %> Integrator coefficient
    Ki = (BnTs / (Xi + 1 / (4 * Xi))) ^ 2 / (Kd * K0);


% TASK
% pay attation to frame structure
f = figure;
SNRSNR = -5 : 5 : 30;
freq_off = -0.4 : 0.01 : 0.4;
for i = 1 : length(SNRSNR)
    RMSE = zeros(1, length(freq_off));
    for k = 1 : length(freq_off)
        Channel_IQ = awgn(IQ_TX_Frame, SNRSNR(i), 'measured');
        Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * freq_off(k));
        [RX_IQ_DM, DM_estimate, DM_filtered] = DM(Channel_IQ, Kp, Ki, IQ_SOF);
        RMSE(k) = norm(DM_filtered(901 : end) - freq_off(k)) / sqrt(100);
    end
    semilogy(freq_off, RMSE, 'DisplayName', string(SNRSNR(i)))
    hold on
end
hold off
legend
xlabel('freq offset')
ylabel('RMSE')
grid on
saveas(f, "1.fig")
%% Receiver with frequency estimator based on Delay and Multiply with D=2
% in feedback scheme
 % normalised frequency
f = figure;
XiXi = 0 : 0.5 : 2;
freq_off = 0.2;
SNR = 30;
for i = 1 : length(XiXi) 
% Configurate the Loop Filter
% =========================================================================
% Loop filter preset
% -------------------------------------------------------------------------
    Xi = XiXi(i);           % detector gainDampingFactor
    BnTs = 2.2;       % Normalized loop bandwidth (Ts = 1 for this problem)
    Kd = 2*pi;          % Phase (not change)
    K0 = 1;             % not change
% =========================================================================
%> Loop filter coefficient calculation
% -------------------------------------------------------------------------
    %todo прописать расчет коэффициентов фильтра
    %> Proportional coefficient
    Kp = 2 * Xi * (BnTs / (Xi + 1 / (4 * Xi))) / (Kd * K0);
    %> Integrator coefficient
    Ki = (BnTs / (Xi + 1 / (4 * Xi))) ^ 2 / (Kd * K0);

% TASK
% pay attation to frame structure
    Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
    Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * freq_off);
    [RX_IQ_DM, DM_estimate, DM_filtered] = DM(Channel_IQ, Kp, Ki, IQ_SOF);
    plot(1 : 50, DM_filtered(1:50), 'DisplayName', string(XiXi(i)))
    hold on
end
hold off
legend
title('variation of Xi')
xlabel('Frame number')
ylabel('DM filtered')
grid on
saveas(f, "2.fig")
%%
f = figure;
XiXi = 0 : 0.5 : 2;
freq_off = -0.4 : 0.01 : 0.4;
SNR = 30;
for i = 1 : length(XiXi) 
% Configurate the Loop Filter
% =========================================================================
% Loop filter preset
% -------------------------------------------------------------------------
    Xi = XiXi(i);           % detector gainDampingFactor
    BnTs = 0.5;       % Normalized loop bandwidth (Ts = 1 for this problem)
    Kd = 2*pi;          % Phase (not change)
    K0 = 1;             % not change
% =========================================================================
%> Loop filter coefficient calculation
% -------------------------------------------------------------------------
    %todo прописать расчет коэффициентов фильтра
    %> Proportional coefficient
    Kp = 2 * Xi * (BnTs / (Xi + 1 / (4 * Xi))) / (Kd * K0);
    %> Integrator coefficient
    Ki = (BnTs / (Xi + 1 / (4 * Xi))) ^ 2 / (Kd * K0);

% TASK
% pay attation to frame structure
    RMSE = zeros(1, length(freq_off));
    for k = 1 : length(freq_off)
        Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
        Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * freq_off(k));
        [RX_IQ_DM, DM_estimate, DM_filtered] = DM(Channel_IQ, Kp, Ki, IQ_SOF);
        RMSE(k) = norm(DM_filtered(901 : end) - freq_off(k)) / sqrt(100);
    end
    semilogy(freq_off, RMSE, 'DisplayName', string(XiXi(i)))
    hold on
end
hold off
legend
xlabel('freq offset')
ylabel('RMSE')
grid on
saveas(f, "21.fig")
%% Receiver with frequency estimator based on Delay and Multiply with D=2
% in feedback scheme
 % normalised frequency
BnTsBnTs = 1 : 0.25 : 2.5;
freq_off = 0.2;
SNR = 30;
f = figure;
for i = 1 : length(BnTsBnTs) 
% Configurate the Loop Filter
% =========================================================================
% Loop filter preset
% -------------------------------------------------------------------------
    Xi = 0.5;           % detector gainDampingFactor
    BnTs = BnTsBnTs(i);       % Normalized loop bandwidth (Ts = 1 for this problem)
    Kd = 2*pi;          % Phase (not change)
    K0 = 1;             % not change
% =========================================================================
%> Loop filter coefficient calculation
% -------------------------------------------------------------------------
    %todo прописать расчет коэффициентов фильтра
    %> Proportional coefficient
    Kp = 2 * Xi * (BnTs / (Xi + 1 / (4 * Xi))) / (Kd * K0);
    %> Integrator coefficient
    Ki = (BnTs / (Xi + 1 / (4 * Xi))) ^ 2 / (Kd * K0);


% TASK
% pay attation to frame structure
    Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
    Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * freq_off);
    [RX_IQ_DM, DM_estimate, DM_filtered] = DM(Channel_IQ, Kp, Ki, IQ_SOF);
    plot(1 : 50, DM_filtered(1:50), 'DisplayName', string(BnTsBnTs(i)))
    hold on
end
hold off
legend
xlabel('Frame number')
title('variation on BnTs')
ylabel('DM filtered')
grid on
saveas(f, "3.fig")
%%
BnTsBnTs = 0.5 : 0.1 : 1.5;
freq_off = -0.4 : 0.01 : 0.4;
SNR = 30;
f = figure;
for i = 1 : length(BnTsBnTs) 
% Configurate the Loop Filter
% =========================================================================
% Loop filter preset
% -------------------------------------------------------------------------
    Xi = 1;           % detector gainDampingFactor
    BnTs = BnTsBnTs(i);       % Normalized loop bandwidth (Ts = 1 for this problem)
    Kd = 2*pi;          % Phase (not change)
    K0 = 1;             % not change
% =========================================================================
%> Loop filter coefficient calculation
% -------------------------------------------------------------------------
    %todo прописать расчет коэффициентов фильтра
    %> Proportional coefficient
    Kp = 2 * Xi * (BnTs / (Xi + 1 / (4 * Xi))) / (Kd * K0);
    %> Integrator coefficient
    Ki = (BnTs / (Xi + 1 / (4 * Xi))) ^ 2 / (Kd * K0);


% TASK
% pay attation to frame structure
    RMSE = zeros(1, length(freq_off));
    for k = 1 : length(freq_off)
        Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
        Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * freq_off(k));
        [RX_IQ_DM, DM_estimate, DM_filtered] = DM(Channel_IQ, Kp, Ki, IQ_SOF);
        RMSE(k) = norm(DM_filtered(901 : end) - freq_off(k)) / sqrt(100);
    end
    semilogy(freq_off, RMSE, 'DisplayName', string(BnTsBnTs(i)))
    hold on
end
hold off
legend
xlabel('freq offset')
ylabel('RMSE')
grid on
saveas(f, "31.fig")
%% Receiver with frequency estimator based on Luise and Reggiannini for N = 20 
% in feedforward scheme
f = figure;
freq_off = -0.4 : 0.01 : 0.4;
SNRSNR = -5 : 5 : 30;
% TASK
% pay attation to frame structure
RMSE = zeros(1, length(freq_off));
for i = 1 : length(SNRSNR)
    RMSE = zeros(1, length(freq_off));
    for k = 1 : length(freq_off)
        Channel_IQ = awgn(IQ_TX_Frame, SNRSNR(i), 'measured');
        Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * freq_off(k));
        [RX_IQ_LR, LR_estimate] = LR(Channel_IQ, IQ_SOF);
        RMSE(k) = norm(LR_estimate(901 : end) - freq_off(k)) / sqrt(100);
    end
    semilogy(freq_off, RMSE, 'DisplayName', string(SNRSNR(i)))
    hold on
end
hold off
legend
xlabel('freq offset')
ylabel('RMSE')
grid on
saveas(f, "4.fig")
%%
f = figure;
Freq_Offset = 0.02; % normalised frequency
SNR = 30; % dB
% TASK
% pay attation to frame structure
Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * Freq_Offset);
[RX_IQ_LR, LR_estimate] = LR(Channel_IQ, IQ_SOF);
semilogy(1 : length(LR_estimate), LR_estimate)
legend
xlabel('Frame number')
ylabel('LR estimate')
grid on
saveas(f, "5.fig")
%%
Freq_Offset = 0.027; % normalised frequency
SNR = 30; % dB
Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
Channel_IQ = Channel_IQ .* exp(1j * 2 * (1 : size(Channel_IQ, 2)) * pi * Freq_Offset);
[RX_IQ_LR, LR_estimate] = LR(Channel_IQ, IQ_SOF);
f = figure();
scatter(real(Channel_IQ(1000, 1 : end)), imag(Channel_IQ(1000, 1 : end)), "filled")
title("Plot before")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, ~] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
saveas(f, "Const_LR_before.fig")
f = figure();
scatter(real(RX_IQ_LR(1000, 1 : end)), imag(RX_IQ_LR(1000, 1 : end)), "filled")
title("Plot after")
xlabel("In-Phase")
ylabel("Quadrature")
xlim([-2 2])
ylim([-2 2])
[Dictionary, ~] = constellation_func("QPSK");
Dict = name("QPSK");
for i = 1:length(Dictionary)
    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
end
axis equal
grid on
saveas(f, "Const_LR_after.fig")
%% Analysis
% =========================================================================
% TASK
% Compare results and make a conclusion
% Compare the RMSE for different frequency offsets for both synchronisation schemes
% Compare the time of synchronisation
% -------------------------------------------------------------------------

