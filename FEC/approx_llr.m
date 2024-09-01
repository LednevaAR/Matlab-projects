function [llr] = approx_llr(IQ, Constellation, snr)
    [Dictionary, D, Bit_depth_Dict] = constellation_func(Constellation);
    D = D';
    llr = zeros(1, length(IQ) * Bit_depth_Dict);
    sigma = 1 / (10 ^ (snr / 10));
    for i = 1 : length(IQ)
        A = (abs(IQ(i) - Dictionary).^2) ./ sigma^2;
        for j = 1 : Bit_depth_Dict
            B = -A .* (D(j, :) - 1);
            C = A .* D(j, :);
            llr((i - 1) * Bit_depth_Dict + j) = -(min(B(B ~= 0)) - min(C(C ~= 0)));
        end
    end
    llr = llr';
end

