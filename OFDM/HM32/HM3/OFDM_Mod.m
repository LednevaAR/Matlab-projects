function Tx_OFDM_symbols = OFDM_Mod(Tx_IQ_points, N_fft, N_carrier, T_guard)
    Tx_OFDM_symbols = zeros(numel(Tx_IQ_points) / N_carrier, N_fft + T_guard);
    for i = 1 : numel(Tx_IQ_points) / N_carrier
        Tx_IQ = Tx_IQ_points((i - 1) * N_carrier + 1 : i * N_carrier);
        Tx_OFDM_symbols(i, T_guard + 1 : end) = ifft(Tx_IQ, N_fft);
    end
    Tx_OFDM_symbols(1 : end, 1 : T_guard) = Tx_OFDM_symbols(1 : end, end - T_guard + 1 : end);
end