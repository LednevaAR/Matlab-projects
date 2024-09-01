function Input_Bit_Buffer = file_reader(File, len)
    Input = imread(File);
    Input_Buffer = reshape(Input, 1, numel(Input));
    Input_Bit_Buffer = dec2bin(Input_Buffer, 8);
    Input_Bit_Buffer = double(reshape(Input_Bit_Buffer, 1, numel(Input_Bit_Buffer))) - 48;
    Input_Bit_Buffer = Input_Bit_Buffer(1 : len);
end

