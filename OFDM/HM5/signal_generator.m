function Tx_OFDM_Signal = signal_generator(Tx_OFDM_symbols)
    Tx_OFDM_Signal = conj(reshape(Tx_OFDM_symbols', 1, numel(Tx_OFDM_symbols)));
end

