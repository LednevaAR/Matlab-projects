function PN = RSLOS(Reg)
    PN = zeros(1, 1e5);
    for i = 1:1e5
        Reg = circshift(Reg, 1);
        Reg(1) = xor(Reg(1), Reg(numel(Reg)));
        PN(i) = Reg(1);
    end
end