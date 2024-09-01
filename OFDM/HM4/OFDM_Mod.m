function Tx_OFDM_symbols = OFDM_Mod(Tx_IQ_points, N_fft, T_guard, inform_index, pilot_index, amp_pilots, phase_mod)
    Tx_OFDM_symbols = zeros(numel(Tx_IQ_points) / length(inform_index), N_fft + T_guard);
    phase_pilots = zeros(1, length(pilot_index));
    if (phase_mod == 1)
        phase_pilots(1 : end) = exp(1j * 2 * pi * 1 / 2);
    elseif (phase_mod == 2)
        phase_pilots(1 : 2 : end) = exp(1j * 2 * pi * 0);
        phase_pilots(2 : 2 : end) = exp(1j * 2 * pi * 1 / 2);
    else
        phase_pilots(1 : end) = exp(1j * 2 * pi * randi([0, 1], 1, length(phase_pilots)) / 2);
    end
    Tx_IQ = zeros(1, length(inform_index) + length(pilot_index));
    for i = 1 : numel(Tx_IQ_points) / length(inform_index)
        Tx_IQ(inform_index) = Tx_IQ_points((i - 1) * length(inform_index) + 1 : i * length(inform_index));
        Tx_IQ(pilot_index) = amp_pilots * max(abs(Tx_IQ_points)) .* phase_pilots;
        Tx_OFDM_symbols(i, T_guard + 1 : end) = ifft(Tx_IQ, N_fft);
    end
    Tx_OFDM_symbols(1 : end, 1 : T_guard) = Tx_OFDM_symbols(1 : end, end - T_guard + 1 : end);
end