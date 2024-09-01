function m_seq = RSLOS(LH, coeff)
    bin = oct2poly(coeff);
    m_seq = zeros(1, LH);
    Reg = zeros(1, length(bin) - 1); % Reg = randi([0 1], length(bin) - 1, 1, 'boolean')';
    for i = 1 : LH
        Reg = circshift(Reg, 1);
        Reg(1) = xor(Reg(1), Reg(numel(Reg)));
        m_seq(i) = Reg(1);
    end
end