function nsignal = Noise (signal, SNR)
 
    sigma = sqrt (1/10^(SNR/10));
    sizeSig = size(signal);

    M = normrnd(0, sigma,[2, sizeSig(1,2)]);
    M = M ./sqrt(2);
    nsignal = zeros (1, sizeSig(1,2));
    for j = 1 : sizeSig(1,2)
        nsignal(j) = signal(j) +  M(1, j)*1i + M(2, j);
    end

end


function NoisedSignal = NoiseGenerator(SNR, Signal)
    Noise = rms(Signal)/(10 ^ (SNR/20)) * normrnd(0, 1, 1, numel(Signal));
    iNoise = 1i * rms(Signal)/(10 ^ (SNR/20)) * normrnd(0, 1, 1, numel(Signal)); 
    NoisedSignal = Signal + Noise + iNoise;
end
