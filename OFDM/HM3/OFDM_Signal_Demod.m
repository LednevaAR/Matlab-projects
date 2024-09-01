function Rx_OFDM_data = OFDM_Signal_Demod(Rx_OFDM_Signal, T_guard, N_fft)
    k = 1;
    for i = 1 : T_guard + N_fft : numel(Rx_OFDM_Signal)
        Rx_OFDM_data(k, 1 : N_fft) = Rx_OFDM_Signal(i + T_guard : i + T_guard + N_fft - 1);
        k = k + 1;
    end
end