function rand = random(Input_Bit_Buffer, Register, Amount_OFDM_Frames)
    rand = zeros(1, numel(Input_Bit_Buffer));
    for i = 1 : numel(Input_Bit_Buffer)
        if mod(i, numel(Input_Bit_Buffer) / Amount_OFDM_Frames) == 1
            Reg = Register;
        end
        Reg = circshift(Reg, 1);
        Reg(1) = xor(Reg(1), Reg(numel(Reg)));
        rand(i) = xor(Reg(1), Input_Bit_Buffer(i));
    end
end