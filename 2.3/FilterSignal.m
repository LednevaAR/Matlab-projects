function FilteredNoisedSignal = FilterSignal(NoisedSignal)
    F = zeros(1, numel(NoisedSignal));
    if (numel(F) >= 42)
        F(7:42) = ones(1, 42 - 6);
    elseif (numel(F) >= 7)
        F(7:end) = ones(1, numel(F) - 6);
    end
    FilteredNoisedSignal = ifft(fft(NoisedSignal) .* F);
end