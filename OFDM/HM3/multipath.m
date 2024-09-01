function result_signal = multipath(Rx_OFDM_Signal, Channel)
    result_signal = zeros(size(Rx_OFDM_Signal));
    for i = 1 : size(Channel, 1)
        result_signal = result_signal + Channel(i, 2) * [Rx_OFDM_Signal(1, Channel(i, 1) + 1 : end), zeros(1, Channel(i, 1))];
    end
end

