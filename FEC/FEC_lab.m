clear; clc; close all;
%%
%============================= Part 1  Warmup =====================================
% configuring LDPC encoder and decoder

% prototype matrix as defnied in Wi-Fi (IEEE® 802.11)
P = [
    16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
    25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
    25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
     9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
    24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
     2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
blockSize = 27; % N.B. this blocksize is connected with optimized pairty check matrix generation method 
                % it is NOT the blocksize of the ldpc
                % https://prezi.com/aqckvai6jux-/ldpc/?utm_campaign=share&utm_medium=copy
H = ldpcQuasiCyclicMatrix(blockSize, P); % getting parity-check matrix

cfgLDPCEnc = ldpcEncoderConfig(H); % configuring encoder
cfgLDPCDec = ldpcDecoderConfig(H); % configuring decoder

% using cfgLDPCEnc variable, print our the number of inofrmation, parity
% check bits and the coderate
fprintf('Number of information bits in a block: %d\n', cfgLDPCEnc.NumInformationBits);
fprintf('Number of parity check bits in a block: %d\n', cfgLDPCEnc.NumParityCheckBits);
coderate = cfgLDPCEnc.NumInformationBits / cfgLDPCEnc.BlockLength;
fprintf('Coderate: %f\n', coderate);

%% simple test to check that encoder and decoder configured correctly

test_message = boolean(randi([0 1],cfgLDPCEnc.NumInformationBits, 1, 'int8'));
encodedData = ldpcEncode(test_message, cfgLDPCEnc);

% calculate the syndrome
s = H * encodedData; %YOUR CODE HERE
s = mod(s, 2); % we need xor instead of multiplication
if(~any(s))
    fprintf('No errors!\n');
else
    fprintf('Errors detected during syndrome check!\n');
end

% deliberately distorting one bit of the message
encodedData(randi(numel(encodedData))) = ~(encodedData(randi(numel(encodedData))));

% checking the syndrome once again
s = H * encodedData; %YOUR CODE HERE
s = mod(s, 2); % we need xor instead of multiplication
if(~any(s))
    fprintf('No errors!\n');
else
    fprintf('Errors detected during syndrome check!\n');
end
%% ============= Part 2 comparing coded and uncoded system =================

maxnumiter = 10;
snr = 0 : 10; % adjust the snr range to the constellation you choose! QPSK
numframes = 10000;

% check manual on the build-in ber counter
% it outputs three variables
ber = comm.ErrorRate; %build-in BER counter
ber2 = comm.ErrorRate; %build-in BER counter

% arrays to store error statistic
errStats = zeros(length(snr), numframes, 3); 
errStatsNoCoding = zeros(length(snr), numframes, 3);
tStart = tic;
for ii = 1:length(snr)
    for counter = 1:numframes
        data = randi([0 1], cfgLDPCEnc.NumInformationBits, 1, 'int8');
        % Transmit and receive with LDPC coding
        encodedData = ldpcEncode(data, cfgLDPCEnc);
        
        % YOUR MAPPER HERE choose any constellation type you like
        modSignal = mapping(encodedData, "QPSK"); %YOUR CODE HERE Mapper(encodedData, ...);

        [rxsig, noisevar] = awgnoise(modSignal, snr(ii)); % use yiur AWGN function

        % YOUR DEMAPPER HERE N.B. Perform Soft Demapping, output llr!
        llr = softdemapping(rxsig, "QPSK", snr(ii));% YOUR CODE HERE Demapper(rxsig, ...);

        rxbits = ldpcDecode(llr, cfgLDPCDec, maxnumiter);
        errStats(ii, counter, :) = ber(data, rxbits);
        %========================================
        
        % no coding system
        noCoding = mapping(data, "QPSK");%YOUR CODE HERE Mapper(encodedData, ...);
        [rxNoCoding, noisevar] = awgnoise(noCoding,snr(ii));
        % YOUR DEMAPPER HERE N.B. Perform Hard Demapping, output bits!
        rxBitsNoCoding = harddemapping(rxNoCoding, "QPSK");% YOUR CODE HERE Demapper(rxNoCoding, ...);
        errStatsNoCoding(ii, counter, :) = ber2(data,int8(rxBitsNoCoding));
    end
    fprintf(['SNR = %2d\n   Coded: Error rate = %1.2f, ' ...
        'Number of errors = %d\n'], ...
        snr(ii),mean(errStats(ii, :, 1), 2), mean(errStats(ii, :, 2), 2))
    fprintf(['Noncoded: Error rate = %1.2f, ' ...
        'Number of errors = %d\n'], ...
        mean(errStatsNoCoding(ii, :, 1), 2), mean(errStatsNoCoding(ii, :, 2), 2))
    reset(ber);
    reset(ber2);
end
ber.release();
ber2.release();
tend = toc(tStart);
fprintf('Simulation finished after %.2f s\n', tend);

%%
f = figure;
semilogy(snr, mean(errStatsNoCoding(:, :, 1), 2), 'LineWidth', 2, 'DisplayName','NoCoding')
hold on
semilogy(snr, mean(errStats(:, :, 1), 2), 'LineWidth', 2, 'DisplayName','LDPC coding')
hold off
legend
xlabel('SNR, dB');
ylabel('BER')
grid on
set(gca, 'Fontsize', 20)
saveas(f, "Part_2_SNR.fig")
% save('BER_SNR_results.mat', 'errStatsNoCoding', 'errStats', '-v7.3')

% Replot the results in BER vs Eb/N0 scale
f = figure;
Eb_N0 = Eb_N0_convert(snr, "QPSK");
semilogy(Eb_N0, mean(errStatsNoCoding(:, :, 1), 2), 'LineWidth', 2, 'DisplayName','NoCoding')
hold on
semilogy(Eb_N0, mean(errStats(:, :, 1), 2), 'LineWidth', 2, 'DisplayName','LDPC coding')
hold off
legend
xlabel('Eb/N0, dB');
ylabel('BER')
grid on
set(gca, 'Fontsize', 20)
saveas(f, "Part_2_Eb_N0.fig")
% how the shape of curves has changed?
% what is the gain in dB?
%%
% +20 points: compare results with llr and approximate llr formulas
maxnumiter = 10;
snr = 0 : 10;
numframes = 10000;
ber = comm.ErrorRate; %build-in BER counter
ber2 = comm.ErrorRate; %build-in BER counter
errStats = zeros(length(snr), numframes, 3); 
errStatsAppr = zeros(length(snr), numframes, 3);
for ii = 1:length(snr)
    for counter = 1:numframes
        data = randi([0 1], cfgLDPCEnc.NumInformationBits, 1, 'int8');
        % Transmit and receive with LDPC coding
        encodedData = ldpcEncode(data, cfgLDPCEnc);
        
        % YOUR MAPPER HERE choose any constellation type you like
        modSignal = mapping(encodedData, "QPSK"); %YOUR CODE HERE Mapper(encodedData, ...);

        [rxsig, noisevar] = awgnoise(modSignal, snr(ii)); % use yiur AWGN function

        % YOUR DEMAPPER HERE N.B. Perform Soft Demapping, output llr!
        llr = softdemapping(rxsig, "QPSK", snr(ii));% YOUR CODE HERE Demapper(rxsig, ...);
        llrappr = approx_llr(rxsig, "QPSK", snr(ii));
        rxbits = ldpcDecode(llr, cfgLDPCDec, maxnumiter);
        rxbitsappr = ldpcDecode(llrappr, cfgLDPCDec, maxnumiter);
        errStats(ii, counter, :) = ber(data, rxbits);
        errStatsAppr(ii, counter, :) = ber2(data, rxbitsappr);
    end
    reset(ber);
    reset(ber2);
end
ber.release();
ber2.release();
%%
f = figure;
semilogy(snr, mean(errStats(:, :, 1), 2), 'LineWidth', 2, 'DisplayName','llr')
hold on
semilogy(snr, mean(errStatsAppr(:, :, 1), 2), 'LineWidth', 2, 'DisplayName','Apprllr')
hold off
legend
xlabel('SNR, dB');
ylabel('BER')
grid on
set(gca, 'Fontsize', 20)
saveas(f, "Part_2_LLR.fig")
%% ================ Part 3: default LDPC with different numbers of iterations =========================

% change the snr range to capture behaviour of coded curves only
snr2 = 4:0.2:6;

maxnumiters = [5, 20]; % we will plot curves for two values of decoding iterations
numframes = 10000;
errStats_it_num = zeros(length(snr2), numframes, 3, numel(maxnumiters));

tStart = tic;

% +10 points for using parfor here and calculating speedup
for ii = 1:length(snr2)
    for m = 1:numel(maxnumiters)
        maxnumiter = maxnumiters(m);
        parfor counter = 1:numframes
            data = randi([0 1],cfgLDPCEnc.NumInformationBits,1,'int8');
            % Transmit and receive with LDPC coding
            encodedData = ldpcEncode(data,cfgLDPCEnc);

            % YOUR MAPPER HERE choose any constellation type you like
            modSignal = mapping(encodedData, "QPSK");%YOUR CODE HERE Mapper(encodedData, ...);

            [rxsig, noisevar] = awgnoise(modSignal,snr2(ii)); % use yiur AWGN function

            % YOUR DEMAPPER HERE N.B. Perform Soft Demapping, output llr!
            llr = softdemapping(rxsig, "QPSK", snr(ii));% YOUR CODE HERE Demapper(rxsig, ...);

            rxbits = ldpcDecode(llr,cfgLDPCDec,maxnumiter);
            br = ber(data,rxbits);
            errStats_it_num(ii, counter, :, m) = br;
        end
        fprintf(['SNR = %2d\n   Coded with %d iterations: Error rate = %1.5f, ' ...
            'Number of errors = %d\n'], ...
            snr2(ii), maxnumiter, mean(errStats_it_num(ii, :, 1, m), 2), mean(errStats_it_num(ii, :, 2, m), 2))
        reset(ber);
    end
end
ber.release();
tend = toc(tStart);
fprintf('Simulation finished in %.2f s\n', tend);
fprintf('With parfor:  s\n')
fprintf('Without parfor: 174.07, 161.52, 406.41, 342.93 s\n')
fprintf('Speedup = \n')

%%
f = figure;
semilogy(snr, mean(errStatsNoCoding(:, :, 1), 2), 'LineWidth', 2)
hold on
for m = 1:numel(maxnumiters)
    semilogy(snr2, mean(errStats_it_num(:, :, 1, m), 2), 'LineWidth', 2)
    hold on
end
hold off
grid on
xlabel('SNR, dB')
ylabel('BER')
xlim([4 6])
set(gca, 'FontSize', 20)
legend('No coding', strcat('LDPC 3/4 ', 32, num2str(maxnumiters(1)), 32,'iterations'), strcat('LDPC 3/4 ', 32, num2str(maxnumiters(2)), 32, 'iterations'));

%save('BER_SNR_results.mat', 'errStats_it_num', '-append')
saveas(f, "Part_3_SNR.fig")
% change the plot to Eb/N0 scale!
f = figure;
Eb_N0 = Eb_N0_convert(snr, "QPSK");
Eb_N02 = Eb_N0_convert(snr2, "QPSK");
semilogy(Eb_N0, mean(errStatsNoCoding(:, :, 1), 2), 'LineWidth', 2)
hold on
for m = 1:numel(maxnumiters)
    semilogy(Eb_N02, mean(errStats_it_num(:, :, 1, m), 2), 'LineWidth', 2)
    hold on
end
hold off
grid on
xlabel('Eb/N0, dB')
ylabel('BER')
xlim([1 3])
set(gca, 'FontSize', 20)
legend('No coding', strcat('LDPC 3/4 ', 32, num2str(maxnumiters(1)), 32,'iterations'), strcat('LDPC 3/4 ', 32, num2str(maxnumiters(2)), 32, 'iterations'));

saveas(f, "Part_3_Eb_N0.fig")

%% ========================= Part 4: diffrent decoding methods with the same max number of iterations

% https://www.mathworks.com/help/comm/ref/ldpcdecode.html
cfgLDPCDec2 = ldpcDecoderConfig(H, 'norm-min-sum'); % configuring second decoder
maxnumiter = 10;
snr2 = 4:0.2:6;
numframes = 10000;

errStats_minsum = zeros(length(snr2), numframes, 3);
errStats_bp = zeros(length(snr2), numframes, 3);

ber = comm.ErrorRate;
ber2 = comm.ErrorRate;

MinSumScalingFactor = 0.7; % task: find the best parameter 
t_min_sum = zeros(length(snr2), 1);
t_bp = zeros(length(snr2), 1);
%%
tStart = tic;
for ii = 1:length(snr2)
    start = tic;
    for counter = 1:numframes
        data = randi([0 1],cfgLDPCEnc.NumInformationBits,1,'int8');
        % Decode with belief propagation
        encodedData = ldpcEncode(data,cfgLDPCEnc);

        % YOUR MAPPER HERE choose any constellation type you like
        modSignal = mapping(encodedData, "QPSK");%YOUR CODE HERE Mapper(encodedData, ...);

        [rxsig, noisevar] = awgnoise(modSignal,snr2(ii)); % use yiur AWGN function
        
        % YOUR DEMAPPER HERE N.B. Perform Soft Demapping, output llr!
        llr = softdemapping(rxsig, "QPSK", snr(ii));% YOUR CODE HERE Demapper(rxsig, ...);
        
        % decode with MinSum
        %rxbits = ldpcDecode(llr,cfgLDPCDec2,maxnumiter, 'MinSumScalingFactor', MinSumScalingFactor);

        %errStats_minsum(ii, counter, :) = ber(data,rxbits);

        % ================================
        % Decode with layered belief propagation

        rxbits = ldpcDecode(llr,cfgLDPCDec,maxnumiter);
        errStats_bp(ii, counter, :) = ber2(data,rxbits);
    end
    t_bp(ii) = toc(start);
    fprintf(['SNR = %2d\n   Min Sum decoding: Error rate = %e, ' ...
        'Number of errors = %d, average time %.4f s\n'], ...
        snr2(ii),mean(errStats_minsum(ii, :, 1), 2), mean(errStats_minsum(ii, :, 2), 2),  t_min_sum(ii))
    fprintf(['BP decoding: Error rate = %e, ' ...
        'Number of errors = %d, average time %.4f s\n'], ...
        mean(errStats_bp(ii, :, 1), 2), mean(errStats_bp(ii, :, 2), 2), t_bp(ii))
    reset(ber);
    reset(ber2);
end
t = toc(tStart);
fprintf('Simulation finished after %.2f s\n', t);

%%
f = figure();
semilogy(snr2, mean(errStats_minsum(:, :, 1), 2))
hold on
semilogy(snr2, mean(errStats_bp(:, :, 1), 2), '--')
hold off

legend('MinSum', 'Belief Propagation')
xlabel('SNR, dB')
ylabel('BER')
grid on
set(gca, 'FontSize', 20)

%save('BER_SNR_results.mat', 'errStats_minsum', 'errStats_bp', '-append')
saveas(f, "Part_4_SNR.fig")
%% Part four: compare the speed of the algorithms

% compare the speed of Belief Propagation and MinSum decoders
f = figure();
semilogy(snr2, t_min_sum, 'DisplayName', 'MinSum')
hold on
semilogy(snr2, t_bp, 'DisplayName', 'Belief Propagation')
hold off

legend
xlabel('SNR, dB')
ylabel('Time, s')
grid on
set(gca, 'FontSize', 20)

%save('BER_SNR_results.mat', 'errStats_minsum', 'errStats_bp', '-append')
saveas(f, "Part_4_time.fig")
%%
%tStart = tic;
cfgLDPCDec2 = ldpcDecoderConfig(H, 'norm-min-sum'); % configuring second decoder
maxnumiter = 10;
snr2 = 4:0.2:6;
numframes = 10000;

ber = comm.ErrorRate;
ber2 = comm.ErrorRate; 
t_min_sum = zeros(length(snr2), length(MinSumScalingFactor));
MinSumScalingFactor = 0.1:0.1:1.0;
errStats_minsum = zeros(length(snr2), numframes, 3, length(MinSumScalingFactor));
for iii = 1:length(MinSumScalingFactor)
    for ii = 1:length(snr2)
        start = tic;
        for counter = 1:numframes
            data = randi([0 1],cfgLDPCEnc.NumInformationBits,1,'int8');
            % Decode with belief propagation
            encodedData = ldpcEncode(data,cfgLDPCEnc);
    
            % YOUR MAPPER HERE choose any constellation type you like
            modSignal = mapping(encodedData, "QPSK");%YOUR CODE HERE Mapper(encodedData, ...);
    
            [rxsig, noisevar] = awgnoise(modSignal,snr2(ii)); % use yiur AWGN function
            
            % YOUR DEMAPPER HERE N.B. Perform Soft Demapping, output llr!
            llr = softdemapping(rxsig, "QPSK", snr(ii));% YOUR CODE HERE Demapper(rxsig, ...);
            
            % decode with MinSum
            rxbits = ldpcDecode(llr,cfgLDPCDec2,maxnumiter, 'MinSumScalingFactor', MinSumScalingFactor(iii));
    
            errStats_minsum(ii, counter, :, iii) = ber(data,rxbits);
    
            % ================================
            % Decode with layered belief propagation
    
            %rxbits = ldpcDecode(llr,cfgLDPCDec,maxnumiter);
            %errStats_bp(ii, counter, :) = ber2(data,rxbits);
        end
    end
    t_min_sum(ii) = toc(start);
    fprintf(['SNR = %2d\n   Min Sum decoding: Error rate = %e, ' ...
        'Number of errors = %d, average time %.4f s\n'], ...
        snr2(ii),mean(errStats_minsum(ii, :, 1), 2), mean(errStats_minsum(ii, :, 2), 2),  t_min_sum(ii))
    fprintf(['BP decoding: Error rate = %e, ' ...
        'Number of errors = %d, average time %.4f s\n'], ...
        mean(errStats_bp(ii, :, 1), 2), mean(errStats_bp(ii, :, 2), 2), t_bp(ii))
    reset(ber);
    reset(ber2);
end
%t = toc(tStart);
fprintf('Simulation finished after %.2f s\n', t);
%%
f = figure();
for m = 1 : length(MinSumScalingFactor)
    semilogy(snr2, mean(errStats_minsum(:, :, 1, m), 2), 'LineWidth', 2, 'DisplayName', string(MinSumScalingFactor(m)))
    hold on
end
hold off

legend
xlabel('SNR, dB')
ylabel('BER')
grid on
set(gca, 'FontSize', 20)

%save('BER_SNR_results.mat', 'errStats_minsum', 'errStats_bp', '-append')
saveas(f, "Part_4_minsum.fig")