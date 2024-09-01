function [H_C, H_L] = f(L, X)
    n = length(X);
    I = (1 : L) + (0 : n - L).';
    B = string(X(I));
    H_L = entropy(B) / L;
    J = (1 : L - 1) + (0 : n - L).';
    %K = L + (0 : n - L).';
    %C = X(K);
    if (L == 2)
        D = X(J);
    else
        D = string(X(J));
    end
    H_C = entropy(B) - entropy(D);
    %H_C = cond_entropy(C, D);
end

