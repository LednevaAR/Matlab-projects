% Фильтрация
%> @file Filtration.m
% =========================================================================
%> @brief Фильтрация
%> @param sign входной сигнал сигнал
%> @param coeff коэффиценты фильтра
%> @param nsamp число выборок на символ
%> @param UpSampFlag [1] -  фильтр с передискретизацией,[0] - фильтр без передискретизации 
%> @return filtsign отфильтрованный сигнал 
% =========================================================================
function filtsign = Filtration(sign, coeff, nsamp, UpSampFlag)
    %> @todo место для вашего кода
    if (UpSampFlag == 1)
        new_sign = zeros(1, length(sign) * nsamp);
        for i = 1 : nsamp : length(new_sign)
            new_sign(i) = sign(idivide(int16(i), int16(nsamp)) + 1);
        end
        filtsign = zeros(1, length(new_sign));
        coeff = [coeff zeros(1, nsamp)];
        for i = 1 : length(new_sign)
            j = idivide(i - 1, int16(nsamp)) * nsamp + 1;
            ind = max(1, j - length(coeff) + nsamp + 1) : j;
            filtsign(i) = filtsign(i) + sum(new_sign(ind) .* coeff(i - ind + 1));
            %for k = 1 : length(coeff) - nsamp
            %    if (i - k + 1 > 0)
            %        filtsign(i : i + nsamp - 1) = filtsign(i : i + nsamp - 1) + new_sign(i - k + 1) * coeff(k : k + nsamp - 1);
            %    end
            %end
        end
    else
        filtsign = zeros(1, length(sign));
        for i = 1 : length(sign)
            ind = max(1, i - length(coeff) + 1) : i;
            filtsign(i) = filtsign(i) + sum(sign(ind) .* coeff(i - ind + 1));
            %for k = 1 : length(coeff)
            %    if (i - k + 1 > 0) 
            %        filtsign(i) = filtsign(i) + sign(i - k + 1) * coeff(k);
            %    end
            %end
        end
    end
end 