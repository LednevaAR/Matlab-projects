function [llr] = softdemapping(IQ, Constellation, snr)
    [Dictionary, D, Bit_depth_Dict] = constellation_func(Constellation);
    D = D';
    llr = zeros(1, length(IQ) * Bit_depth_Dict);
    sigma = 1 / (10 ^ (snr / 10));
    for i = 1 : length(IQ)
        A = exp(-(abs(IQ(i) - Dictionary).^2) ./ sigma^2);
        for j = 1 : Bit_depth_Dict
            llr((i - 1) * Bit_depth_Dict + j) = log(sum(-A .* (D(j, :) - 1)) / sum(A .* D(j, :)));
        end
    end
    llr = llr';
end

