function Probability = Error_check(Input_Bit_Buffer, Output_Bit_Buffer)
    Probability = sum(Input_Bit_Buffer == Output_Bit_Buffer) / numel(Input_Bit_Buffer);
    %disp(Input_Bit_Buffer - Output_Bit_Buffer)
end

