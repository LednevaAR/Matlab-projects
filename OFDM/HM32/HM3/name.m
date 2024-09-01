function Dict = name(Constellation)
    switch Constellation
        case "BPSK"
            Dict = ["0", "1"];
        case "QPSK"
            Dict = ["00", "01", "10", "11"];
        case "8PSK"
            Dict = ["000", "001", "010", "011", "100", "101", "110", "111"];
        case "16-QAM"
            Dict = ["0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"];
    end
end
