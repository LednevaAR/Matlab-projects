function x = MP(y, M, epsilon, max_iter)
    r = y;
    x = zeros(1, size(y, 2));
    r_prev = y;
    I = [];
    M1 = M;
    for n = 1 : max_iter
        [~, idx] = max(abs(r * M1));
        I = [I idx];
        r_prev = r;
        r = r_prev - dot(M1(:, idx), r_prev.') * M1(:, idx);
        x(idx) = dot(M1(:, idx), r_prev.');
        if (sqrt(sum(abs(r - r_prev) .^ 2))) / (sqrt(sum(abs(r_prev) .^ 2))) < epsilon      
            break;
        end
        M1(:, idx) = [];
    end
end