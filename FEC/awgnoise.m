function [rx, noisevar] = awgnoise(signal, snr)
    Noise = 1/(10 ^ (snr/20)) * normrnd(0, 1, 1, numel(signal))/sqrt(2);
    iNoise = 1i * 1/(10 ^ (snr/20)) * normrnd(0, 1, 1, numel(signal))/sqrt(2); 
    rx = signal + Noise + iNoise;
    noisevar = var(Noise + iNoise);
end

