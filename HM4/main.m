clear; clc; close all;
%%
r = 7;
LH = 2^7 - 1;
LD = 10 * LH;
N = 1; % amount of frames
%%
coeff = 217;
m_seq = RSLOS(LH, coeff);
disp(sum(m_seq)) % число единиц
disp(abs(sum(m_seq - 1))) % число нулей
Head = mapping(m_seq, "BPSK");
AutoCorrm = AutoCorr(Head);
BN = linspace(1, numel(AutoCorrm), numel(AutoCorrm));
f = figure();
plot(BN, AutoCorrm);
xlabel('Bit number');
ylabel('AutoCorrelation function');
title('AutoCorrelation function of m-seq');
saveas(f, 'ACF_m_seq.fig')
%%
Data = randi([0 1], LD * N, 1, 'double')';
Data = reshape(Data, N, []);
Header = repmat(m_seq, N, 1);
Bits = Data;
Bits = reshape(Bits, 1, []);
Signal = mapping(Bits, "QPSK");
Head = mapping(Header, "BPSK");
Signal = [1 1 1 1 1 1 Head Signal];
%%
f = figure();
Freq_Offset = -0.2 : 0.1 : 0.2;
SNR = -10 : 1 : 5;
for i = 1 : length(Freq_Offset)
    R = zeros(length(SNR), length(Signal));
    P = zeros(1, length(SNR));
    for j = 1 : length(SNR)
        for k = 1 : 100
            Signal1 = awgn(Signal, SNR(j), 'measured');
            Signal1 = Signal1 .* exp(-1j .* 2 .* (1 : length(Signal1)) .* pi * Freq_Offset(i));
            y = Signal1;
            for n = 1 : length(y)
                sum = 0;
                sum1 = 0;
                sum2 = 0;
                for k = 1 : LH
                    if (n + k - 1 <= length(y))
                        sum = sum + y(n + k - 1) * conj(Head(k));
                        sum1 = sum1 + abs(y(n + k - 1)) ^ 2;
                    end
                    sum2 = sum2 + abs(Head(k)) ^ 2;
                end
                R(j, n) = abs(sum / sqrt(sum1 * sum2));
            end
            [~, index] = max(R(j, 1 : end));
            if (index == 7)
                P(j) = P(j) + 1;
            end
        end
    end
    plot(SNR, P ./ 100, 'LineWidth', 2, 'DisplayName', 'FO=' + string(Freq_Offset(i)) + 'f')
    hold on
end
hold off
legend
grid on
set(gca, 'FontSize', 20)
saveas(f, "Image0.fig")
%%
f = figure(); 
Freq_Offset = -0.2 : 0.1 : 0.2;
SNR = -10 : 1 : 5;
for i = 1 : length(Freq_Offset)
    R = zeros(length(SNR), length(Signal));
    P = zeros(1, length(SNR));
    for j = 1 : length(SNR)
        for k = 1 : 100
            Signal1 = awgn(Signal, SNR(j), 'measured');
            Signal1 = Signal1 .* exp(-1j .* 2 .* (1 : length(Signal1)) .* pi * Freq_Offset(i));
            y = Signal1;
            y = [y 0];
            cx = zeros(1, LH);
            cy = zeros(1, length(y));
            for n = 1 : length(y)
                if (n < length(y))
                    cy(n) = y(n) * conj(y(n + 1));
                end
            end
            for k = 1 : LH
                if (k < LH) 
                    cx(k) = Head(k) * conj(Head(k + 1));
                end
            end
            for n = 1 : length(y)
                sum = 0;
                sum1 = 0;
                sum2 = 0;
                for k = 1 : LH
                    if (n + k - 1 <= length(y))
                        sum = sum + cy(n + k - 1) * conj(cx(k));
                        sum1 = sum1 + abs(cy(n + k - 1)) ^ 2;
                    end
                    sum2 = sum2 + abs(cx(k)) ^ 2;
                end
                R(j, n) = abs(sum / sqrt(sum1 * sum2));
            end
            [~, index] = max(R(j, 1 : end));
            if (index == 7)
                P(j) = P(j) + 1;
            end
        end
    end
    plot(SNR, P ./ 100, 'LineWidth', 2, 'DisplayName', 'FO=' + string(Freq_Offset(i)) + 'f')
    hold on
end
hold off
legend
grid on
set(gca, 'FontSize', 20)
saveas(f, "Image0diff.fig")
%%
%100 кадров, искать максимальный пик, индекс, затем сравнивать с настоящим положением