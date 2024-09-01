function [RX_IQ_LR, LR_E] = LR(Channel_IQ, IQ_SOF)
    LR_NCO = 0;
    N = 20;
    LR_E = zeros(1, size(Channel_IQ, 1));
    RX_IQ_LR = zeros(size(Channel_IQ, 1), size(Channel_IQ, 2));
    for itter_time = 1 : size(Channel_IQ, 1)
        % TASK
        % LR detector
        RX_IQ = Channel_IQ(itter_time, 1 : end);
        z = RX_IQ(1 : length(IQ_SOF)) .* conj(IQ_SOF);
        LR_estimate = (1 / (N * pi)) * angle(SumAutoCorr(z, N)); 
        % TASK
        LR_E(itter_time) = LR_estimate;
        % NCO | Phase Accumulation
        LR_NCO = LR_estimate * (1 : length(RX_IQ));
        % TASK    
        % Compensation
        RX_IQ_LR(itter_time, 1 : end) = RX_IQ .* exp(-1j * 2 * pi * LR_NCO);
    end

% =========================================================================
% How does the estimate behave? Show on the plot
% How did the constellation change?
% -------------------------------------------------------------------------
    %f = figure();
    %scatter(real(RX_IQ_LR), imag(RX_IQ_LR), "filled")
    %title("Plot before")
    %xlabel("In-Phase")
    %ylabel("Quadrature")
    %xlim([-2 2])
    %ylim([-2 2])
    %[Dictionary, Bit_depth_Dict] = constellation_func("QPSK");
   % Dict = name("QPSK");
    %for i = 1:length(Dictionary)
    %    text(real(Dictionary(i)), imag(Dictionary(i)), '\leftarrow' + Dict(i));
    %end
    %axis equal
    %grid on
    %saveas(f, "Const_3_before.fig")
end

