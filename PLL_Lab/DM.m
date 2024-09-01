function [RX_IQ_DM, DM_E, DM_F] = DM(Channel_IQ, Kp, Ki, IQ_SOF)
    D = 2;
    DM_NCO = zeros(1, size(Channel_IQ, 2));
    DM_F = zeros(1, size(Channel_IQ, 1));
    DM_E = zeros(1, size(Channel_IQ, 1));
    p = 0;
    RX_IQ_DM = zeros(size(Channel_IQ, 1), size(Channel_IQ, 2));
    for itter_time = 1 : size(Channel_IQ, 1)
        % Compensation
        RX_IQ = Channel_IQ(itter_time, 1 : end) .* exp(-1j * 2 * pi * DM_NCO);
        % DM detector
        z = RX_IQ(1 : length(IQ_SOF)) .* conj(IQ_SOF);
        DM_estimate = (1 / (2 * pi * D)) * angle(sum(z(D + 1 : end) .* conj(z(1 : end - D))));
        % Loop filter
        DM_filtered = Kp * DM_estimate + Ki * DM_estimate + p;
        % NCO | Phase Accumulation
        DM_NCO = DM_filtered * (1 : length(RX_IQ));
        RX_IQ_DM(itter_time, 1 : end) = RX_IQ;
        DM_F(itter_time) = DM_filtered;
        DM_E(itter_time) = DM_estimate;
        p = DM_filtered - Kp * DM_estimate;
    end
% =========================================================================
% TASK
% For different Damping Factor and BnTs calculate coefficients of loop filter
% What changes in synchronisation when the loop filter coefficients are recalculated?
% Illustrate these changes on the graphs
% How did the constellation change?
% -------------------------------------------------------------------------
   % f = figure();
    %scatter(real(RX_IQ_DM), imag(RX_IQ_DM), "filled")
    %title("Plot before")
    %xlabel("In-Phase")
    %ylabel("Quadrature")
    %xlim([-2 2])
    %ylim([-2 2])
    %[Dictionary, Bit_depth_Dict] = constellation_func("QPSK");
   % Dict = name("QPSK");
   % for i = 1:length(Dictionary)
   %     text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
   % end
   % axis equal
   % grid on
    %saveas(f, "Const_3_before.fig")
end

