% Генерация коэффицентов (импульсной характеристики) для фильтра 
% приподнятого косинуса
%> @file RCcoeff.m
% =========================================================================
%> @brief Генерация коэффицентов (импульсной характеристики)для фильтра 
%> приподнятого косинуса
%> @param span Длина фильтра в символах (количество боковых лепестков sinc, 
%> сумма с двух сторон)
%> @param nsamp Число выборок на символ
%> @param rolloff Коэффицент сглаживания (alfa)
%> @return coeff Коэффиценты для фильтра приподнятого косинуса
% =========================================================================
function coeff = RCcoeff (span, nsamp, rolloff)
    %> @todo Место для вашего кода
    coeff = zeros(1, span * nsamp + 1);
    for i = 1 : size(coeff, 2)
        t = (i - span * nsamp / 2 - 1) / nsamp;
        t_sym = 1;
        if (t == 0)
            coeff(i) = 1 / sqrt(t_sym);
        elseif (abs(t) == t_sym / (2 * rolloff))
            coeff(i) = (1 / sqrt(t_sym)) * (rolloff / 2) * sin(pi / (2 * rolloff));
        else
            coeff(i) = (1 / sqrt(t_sym)) * (sin(pi * t / t_sym) / (pi * t / t_sym)) * cos(pi * t * rolloff / t_sym) / (1 - (2 * rolloff * t / t_sym)^2);
        end
    end
    coeff = coeff / sqrt(sum(coeff .^ 2));
end
