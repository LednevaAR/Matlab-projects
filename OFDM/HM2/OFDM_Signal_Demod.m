function Rx_OFDM_symbols = OFDM_Signal_Demod(Rx_OFDM_Signal, T_guard, N_fft)
    %1ый способ
%     i = 1;
%     k = 1;
%     p = 0;
%     size = 0;
%     while (i <= numel(Rx_OFDM_Signal) - T_guard)
%         if (Rx_OFDM_Signal(i : i + T_guard - 1) == zeros(1, T_guard))
%             if (p ~= 0)
%                 if (k == 1)
%                     size = i - p;
%                 end
%                 Rx_OFDM_symbols(k, 1 : size) = Rx_OFDM_Signal(p : i - 1);
%                 k = k + 1;
%             end
%             p = i + T_guard;
%         end
%         i = i + 1;
%     end
%     Rx_OFDM_symbols(k, 1 : size) = [Rx_OFDM_Signal(p : p + size - 1)];
    %другой способ - можно через N_fft
    k = 1;
    for i = 1 : T_guard + N_fft : numel(Rx_OFDM_Signal)
        Rx_OFDM_symbols(k, 1 : N_fft) = Rx_OFDM_Signal(i + T_guard : i + T_guard + N_fft - 1);
        k = k + 1;
    end
end