function [IQ] = mapping(Bit_Tx, Constellation)
% Make the different dictionary for BPSK, QPSK, 8PSK, 16QAM constellations
% calculate the Bit_depth for each contellation
    [Dictionary, ~, Bit_depth_Dict] = constellation_func(Constellation);
% write  the function of mapping from bit vector to IQ vector
    IQ = zeros(1, ceil(length(Bit_Tx) / Bit_depth_Dict));
    for i = 1:Bit_depth_Dict:length(Bit_Tx) - Bit_depth_Dict + 1
        index = 0;
        for j = 0:Bit_depth_Dict - 1
            index = index + (2^(Bit_depth_Dict - j - 1)) * Bit_Tx(i + j);
        end
        IQ(round((i - 1)/Bit_depth_Dict) + 1) = Dictionary(index + 1);
    end
end



