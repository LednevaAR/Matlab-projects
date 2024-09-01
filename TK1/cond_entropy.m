function H = cond_entropy(X, Y)
    H = 0;
    [A, P] = alphabet_probabilities(Y);
    for i = 1 : length(A)
        H = H + P(i) * cond_val_entropy(X, Y, A(i));
    end
end

