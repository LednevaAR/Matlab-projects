% Генерация коэффицентов (импульсной характеристики)для фильтра корень 
% из приподнятого косинуса
%> @file sqRCcoeff.m
% =========================================================================
%> @brief Генерация коэффицентов (импульсной характеристики)для фильтра корень 
%> из приподнятого косинуса
%> @param span Длина фильтра в символах (количество боковых лепестков sinc, 
%> сумма с двух сторон)
%> @param nsamp Число выборок на символ
%> @param rolloff Коэффицент сглаживания (alfa)
%> @return coeff коэффиценты для фильтра корень из приподнятого косинуса
% =========================================================================
function coeff = sqRCcoeff(span, nsamp, rolloff)
    %> @todo Место для вашего кода
    coeff = zeros(1, span * nsamp + 1);
    for i = 1 : size(coeff, 2)
        %if (abs(f) <= (1 - rolloff) / 2)
        %    coeff(i) = 1;
        %elseif (abs(f) >= (1 + rolloff) / 2)
        %    coeff(i) = 0;
        %else
        %    coeff(i) = sqrt(1 / 2 + (1 / 2) * sin(pi * (1 / 2 - abs(f)) / rolloff));
        %end
        t = (i - span * nsamp / 2 - 1) / nsamp;
        t_sym = 1;
        if (t == 0)
            coeff(i) = (1 / sqrt(t_sym)) * (1 + rolloff * (4 / pi - 1));
        elseif (abs(t) == t_sym / (4 * rolloff))
            coeff(i) = (1 / sqrt(t_sym)) * (rolloff / sqrt(2)) * ((1 + 2 / pi) * sin(pi / (4 * rolloff)) + (1 - 2 / pi) * cos(pi / (4 * rolloff)));
        else
            coeff(i) = (1 / sqrt(t_sym)) * ((4 * rolloff / pi) * cos((1 + rolloff) * pi * t / t_sym) + (t_sym / (pi * t)) * sin((1 - rolloff) * pi * t / t_sym)) / (1 - 16 * rolloff^2 * t^2 / t_sym^2);
        end
    end
    coeff = coeff / sqrt(sum(coeff .^ 2));
end