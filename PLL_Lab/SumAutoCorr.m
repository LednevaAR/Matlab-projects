function SumCorr = SumAutoCorr(PN, N)
    SumCorr = 0;
    for i = 1 : N - 1
        R = sum(PN(i + 1 : end) .* conj(PN(1 : end - i)));
        SumCorr = SumCorr + R / (length(PN) - i);
    end
end