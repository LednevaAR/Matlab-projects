function [Bit] = demapping(IQ, Constellation)
% Make the different dictionary for BPSK, QPSK, 8PSK, 16QAM constellations
% calculate the Bit_depth for each contellation
    [Dictionary, Bit_depth_Dict] = constellation_func(Constellation);
% write  the function of mapping from IQ vector to bit vector
    Bit = zeros(1, length(IQ) * Bit_depth_Dict);
    for i = 1:length(IQ)
        [M, index] = min(abs(Dictionary - IQ(i)));
        index = index - 1;
        for j = Bit_depth_Dict:-1:1
            Bit((i - 1) * Bit_depth_Dict + j) = mod(index, 2);
            index = (index - mod(index, 2))/2;
        end
    end
end


