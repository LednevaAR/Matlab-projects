function H = cond_val_entropy(X, Y, y)
    n = (y == Y);
    H = -1;
    if (~isempty(nonzeros(n)))
        Z = X(n == 1);
        [~, R] = alphabet_probabilities(Z);
        P = (R .* length(Z) ./ length(X)) ./ (length(nonzeros(n)) / length(Y));
        H = -sum(P .* log2(P));
    end
end

