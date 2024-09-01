function [PAPR_slid_wind] = compute_PAPR_slid_wind(Tx_OFDM_Signal, size)
    PAPR_slid_wind = zeros(1, numel(Tx_OFDM_Signal) - size + 1);
    for i = 1 : numel(Tx_OFDM_Signal) - size + 1
        PAPR_slid_wind(i) = 10 * log10(max(abs(Tx_OFDM_Signal(i : i + size - 1)) .^ 2) / mean(abs(Tx_OFDM_Signal(i : i + size - 1)) .^ 2));
    end
end

