function [Eb_N0] = Eb_N0_convert(SNR, Constellation)
    switch Constellation
        case "BPSK"
            Bit_depth_Dict = 1;
        case "QPSK"
            Bit_depth_Dict = 2;
        case "8PSK"
            Bit_depth_Dict = 3;
        case "16-QAM"
            Bit_depth_Dict = 4;
    end
     Eb_N0 = SNR - 10*log10(Bit_depth_Dict);

%      Eb_N0 = SNR + 10*log10(Bit_depth_Dict);
end

