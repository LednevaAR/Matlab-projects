function NoisedSignal = NoiseGenerator(SNR, Signal)
    Noise = 1/(10 ^ (SNR/20)) * normrnd(0, 1, 1, numel(Signal))/sqrt(2) * mean(abs(Signal));
    iNoise = 1i * 1/(10 ^ (SNR/20)) * normrnd(0, 1, 1, numel(Signal))/sqrt(2) * mean(abs(Signal)); 
    NoisedSignal = Signal + Noise + iNoise;
end
