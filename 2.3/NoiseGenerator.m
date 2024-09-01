function NoisedSignal = NoiseGenerator(SNR, Signal)
    Noise = rms(Signal)/(10 ^ (SNR/20)) * normrnd(0, 1, 1, numel(Signal));
    NoisedSignal = Signal + Noise;
end

