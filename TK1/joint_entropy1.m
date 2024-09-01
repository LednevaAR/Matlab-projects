function H = joint_entropy1(X,Y)
    H = 0;
    [Alpha, ~] = alphabet_probabilities(Y);
    for i = 1 : length(Alpha)
        n = (Alpha(i) == Y);
        Z = X(n == 1);
        [~, R] = alphabet_probabilities(Z);
        P = R .* length(Z) ./ length(X);
        H = H - sum(P .* log2(P));
    end          
end

