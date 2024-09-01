function result_signal = multipath(Rx_OFDM_Signal, Channel)
    result_signal = zeros(size(Rx_OFDM_Signal));
    for i = 1 : size(Channel, 1)
        result_signal = result_signal + Channel(i, 2) .* [Rx_OFDM_Signal(1, Channel(i, 1) + 1 : end), zeros(1, Channel(i, 1))];
    end
    N_fft = 1024;
    T_guard = N_fft / 8;
    N_carrier = 400;
    Rx_OFDM_data = OFDM_Signal_Demod(result_signal, T_guard, N_fft);
    Rx_IQ = zeros(size(Rx_OFDM_data, 1), N_fft);
    for i = 1 : size(Rx_OFDM_data, 1)
        Rx_IQ(i, 1 : N_fft) = fft(Rx_OFDM_data(i, 1 : end), N_fft);
        Rx_IQ(i, 1 : N_fft) = [Rx_IQ(i, 1 : N_carrier) zeros(1, N_fft - N_carrier)];
    end
    Rx_IQ_points = conj(reshape(Rx_IQ(:, 1 : N_carrier)', 1, numel(Rx_IQ(:, 1 : N_carrier))));
    Rx_OFDM_symbols = OFDM_Mod(Rx_IQ_points, N_fft, N_carrier, T_guard);
    result_signal = signal_generator(Rx_OFDM_symbols);
end

