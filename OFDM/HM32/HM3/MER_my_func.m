function [MER] = MER_my_func(IQ_RX, Constellation)
    [Dictionary, Bit_depth_Dict] = constellation_func(Constellation);
    sum1 = 0;
    sum2 = 0;
    for i = 1:length(IQ_RX)
        [M, index] = min(abs(Dictionary - IQ_RX(i)));
        sum1 = sum1 + ((real(Dictionary(index)))^2 + (imag(Dictionary(index)))^2);
        sum2 = sum2 + ((real(IQ_RX(i) - Dictionary(index)))^2 + (imag(IQ_RX(i) - Dictionary(index)))^2);
    end
    MER = 10 * log10(sum1/sum2);
end

