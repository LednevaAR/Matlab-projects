function H = entropy(symb_arr)
    [~, emp_prob] = alphabet_probabilities(symb_arr);
    H = -sum(emp_prob .* log2(emp_prob));
end

