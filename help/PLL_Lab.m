%%
close all; clear; clc;

%% config
Freq_Offset = -0.2; % normalised frequency
SNR = 0; % 30 dB 2.5
Eb_N0=Eb_N0_convert(SNR,'QPSK');

%% Transmitter
% config
Amount_of_Frame = 1000;
Length_Data_IQ = 1440;

% Start of frame
SOF = [1 0 0 1 1 1 0 1 0 1 0 1 0 1 1 0 0 1 0 0]; 
IQ_SOF = mapping(SOF, 'BPSK'); % Use this sequence on the Rx as a Pilot-Signal

% bit generator

% QAM | mapper
% Generate vector of bits with length = Amount_of_Frame*Length_Data_IQ*2
Tx_Bits = randi([0 1], Amount_of_Frame * Length_Data_IQ * 2, 1)'; % генерация бит; 
% Use your function Mapping.m from previous courses
TX_IQ_Data = mapping(Tx_Bits, 'QPSK');

% Frame structure 
% |20 IQ BPSK Start-of-Frame| 1440 IQ QPSK|
IQ_TX_Frame = FrameStruct(TX_IQ_Data, IQ_SOF, Amount_of_Frame);



%% Channel
% % Add white Gaussian noise to signal
% Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
% 
% % Add frequency offset
% Channel_IQ=Channel_IQ.*exp(1j.*2.*(1:length(Channel_IQ)).*pi*Freq_Offset);

%% Receiver with frequency estimator based on Delay and Multiply with D=2
% in feedback scheme

% Configurate the Loop Filter
% =========================================================================
% Loop filter preset
% -------------------------------------------------------------------------
    Xi = 1;             %1 detector gainDampingFactor
    BnTs = 0.001;       %0.001 Normalized loop bandwidth (Ts = 1 for this problem)
    Kd = 2*pi;          % Phase (not change)
    K0 = 1;             % not change
% =========================================================================
%> Loop filter coefficient calculation
% -------------------------------------------------------------------------
    %todo прописать расчет коэффициентов фильтра
    wp=BnTs/(Xi+1/(4*Xi));
    %> Proportional coefficient
    Kp = 2*Xi*wp/K0/Kd;
    %> Integrator coefficient
    Ki = wp^2/K0/Kd;


% [RX_IQ_DM, DM_estimate, DM_filtered] = DM(Channel_IQ, Kp, Ki);
% 
% RX_row=RX_IQ_DM(end,:);
% 
% figure();
% scatter(real(RX_row),imag(RX_row));
% figure();
% plot(DM_filtered);
% disp("end");


% TASK
% pay attation to frame structure

Freq_Offset_range = -0.5:0.01:0.5; 

% Initialize arrays to store results
NormRMSFreqEstError_DM = zeros(size(Freq_Offset_range));
NormRMSFreqEstError_LR = zeros(size(Freq_Offset_range));

% Loop through different frequency offsets
SNR_ar=0:5:40;
for j=1:length(SNR_ar)
    SNR=SNR_ar(j);
    for i = 1:length(Freq_Offset_range)
    Freq_Offset = Freq_Offset_range(i);

    % Channel
    Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
    Channel_IQ = Channel_IQ.*exp(1j*2*(1:size(Channel_IQ,2))*pi*Freq_Offset);

    % Receiver with frequency estimator based on DM and LR
    [RX_IQ_DM, DM_estimate, DM_filtered] = DM(Channel_IQ, Kp, Ki);
    [RX_IQ_LR, LR_estimate] = LR(Channel_IQ);

    % Compute normalized RMS frequency estimation error
    NormRMSFreqEstError_DM(i) = norm((DM_filtered(end-500:end)) - Freq_Offset) / sqrt(length((DM_filtered(end-500:end))));
    NormRMSFreqEstError_LR(i) = norm((LR_estimate(end-500:end)) - Freq_Offset) / sqrt(length((LR_estimate(end-500:end))));
    
%     if i==50
%         figure(3);
%         scatter(real(RX_IQ_DM(end,:)),imag(RX_IQ_DM(end,:)))
%         figure(4);
%         scatter(real(RX_IQ_LR(end,:)),imag(RX_IQ_LR(end,:)))
%         disp(1);
%     end
    
    end

    % Plot Norm. RMS Freq. Est. Error vs Norm. Freq. Offset
    figure(1);
    semilogy(Freq_Offset_range, NormRMSFreqEstError_DM, '-o','LineWidth', 2,'MarkerSize', 2);
    xlabel('Normalized Frequency Offset');
    ylabel('Normalized RMS Frequency Estimation Error');
    title('D&M Algorithm SNR 0:5:40');
    grid on;
    hold on;

    figure(2);
    semilogy(Freq_Offset_range, NormRMSFreqEstError_LR, '-s','LineWidth', 2,'MarkerSize', 2);
    xlabel('Normalized Frequency Offset');
    ylabel('Normalized RMS Frequency Estimation Error');
    title('L&R Algorithm for SNR 0:5:40');
    grid on;
    hold on;
end
figure(1);
legend(num2str(SNR_ar'),'Location','southwest')
figure(2);
legend(num2str(SNR_ar'),'Location','southwest')
%% Receiver with frequency estimator based on Luise and Reggiannini for N = 20 
% in feedforward scheme

% TASK
% pay attation to frame structure

% [RX_IQ_LR, LR_estimate] = LR(Channel_IQ);


%% Analysis
% =========================================================================
% TASK
% Compare results and make a conclusion
% Compare the RMSE for different frequency offsets for both synchronisation schemes
% Compare the time of synchronisation
% -------------------------------------------------------------------------

% % Define a range of normalized frequency offsets
% Freq_Offset_range = -0.5:0.1:0.5; 
% 
% % Initialize arrays to store results
% NumSimulations = 10; % adjust as needed
% NormRMSFreqEstError = zeros(length(Freq_Offset_range), NumSimulations);
% NormFreqOffset = zeros(length(Freq_Offset_range), NumSimulations);
% 
% parfor sim = 1:NumSimulations
%     tempNormRMSFreqEstError = zeros(length(Freq_Offset_range), 1);
% 
%     for i = 1:length(Freq_Offset_range)
%         Freq_Offset = Freq_Offset_range(i);
% 
%         % Channel
%         Channel_IQ = awgn(IQ_TX_Frame, SNR, 'measured');
%         Channel_IQ = Channel_IQ .* exp(-1j * 2 * pi * Freq_Offset * (1:size(Channel_IQ,2)));
% 
%         % Receiver with frequency estimator based on Delay and Multiply with D=2
%         [RX_IQ_DM, DM_estimate] = DM(Channel_IQ, Kp, Ki);
% 
%         % Compute normalized RMS frequency estimation error
%         tempNormRMSFreqEstError(i) = norm(DM_estimate - Freq_Offset) / sqrt(length(DM_estimate));
%     end
% 
%     % Assign temporary variables to sliced variables
%     NormRMSFreqEstError(:, sim) = tempNormRMSFreqEstError;
% end
% 
% % Calculate the mean across simulations
% MeanNormRMSFreqEstError = mean(NormRMSFreqEstError, 2);
% 
% % Plot Norm. RMS Freq. Est. Error vs Norm. Freq. Offset
% figure;
% semilogy(Freq_Offset_range, MeanNormRMSFreqEstError, '-o');
% xlabel('Normalized Frequency Offset');
% ylabel('Mean Normalized RMS Frequency Estimation Error');
% title('D&M Algorithm: Mean Norm. RMS Freq. Est. Error vs Norm. Freq. Offset');
% grid on;