function Tx_OFDM_Signal = signal_generator(Tx_OFDM_symbols, T_guard)
    Tx_OFDM_Signal = zeros(size(Tx_OFDM_symbols, 1), size(Tx_OFDM_symbols, 2) + T_guard);
    Tx_OFDM_Signal(1 : end, T_guard + 1 : end) = Tx_OFDM_symbols;
    Tx_OFDM_Signal(1 : end, 1 : T_guard) = Tx_OFDM_symbols(1 : end, end - T_guard + 1 : end);
    Tx_OFDM_Signal = reshape(Tx_OFDM_Signal', 1, numel(Tx_OFDM_Signal));
end

