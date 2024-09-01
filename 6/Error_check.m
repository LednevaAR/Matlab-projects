function [BER] = Error_check(Bit_Tx, Bit_Rx)
    BER = length(find(xor(Bit_Tx(1:length(Bit_Rx)), Bit_Rx) ~= 0))/length(Bit_Rx);
end

