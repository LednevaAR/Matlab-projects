function [Eb_N0] = Eb_N0_convert(SNR, Constellation)
    switch Constellation
        case "BPSK"
            k = 1;
        case "QPSK"
            k = 2;
        case "8PSK"
            k = 3;
        case "16-QAM"
            k = 4;
    end
    Eb_N0 = SNR - 10 * log10(k);
end

