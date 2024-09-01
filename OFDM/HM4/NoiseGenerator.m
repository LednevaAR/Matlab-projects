function NoisedSignal = NoiseGenerator(SNR, Signal)
    Noise = 1/(10 ^ (SNR/20)) * normrnd(0, 1, 1, numel(Signal))/sqrt(2) * sqrt(mean(abs(Signal) .^ 2));
    iNoise = 1i * 1/(10 ^ (SNR/20)) * normrnd(0, 1, 1, numel(Signal))/sqrt(2) * sqrt(mean(abs(Signal) .^ 2)); 
    NoisedSignal = Signal + Noise + iNoise;
end
