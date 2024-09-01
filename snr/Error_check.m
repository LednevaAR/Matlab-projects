function [BER] = Error_check(Bit_Tx, Bit_Rx)

    BER = sum(xor(Bit_Tx,Bit_Rx))/length(Bit_Tx);
end

