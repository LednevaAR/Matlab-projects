function [alphabet, emp_prob] = alphabet_probabilities(symb_arr)
    [alphabet, ~, index] = unique(symb_arr);
    %emp_prob = hist(double(symb_arr), double(alphabet));% ./ length(symb_arr);
    %D = double(symb_arr);
    %emp_prob = (nonzeros(accumarray(D(:), 1, [], @sum)))' ./ length(symb_arr);
    emp_prob = histcounts(index, 1:(length(alphabet) + 1)) ./ length(symb_arr);
end