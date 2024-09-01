function m_seq = RSLOS(LH, coeff)
    bin = oct2poly(coeff);
    disp(bin)
    bin = bin(2 : end);
    m_seq = zeros(1, length(bin) - 1);
    m_seq(length(bin)) = 1;
    % m_seq(1 : length(bin)) = randi([0 1], length(bin), 1, 'double')';
    for i = length(bin) + 1 : length(bin) + LH
        m = 0;
        for j = 1 : length(bin)
            m = xor(m, bin(j) * m_seq(i - j));
        end
        m_seq(i) = m;
    end
    m_seq = m_seq(length(bin) + 1 : end);
end