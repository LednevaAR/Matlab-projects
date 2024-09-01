function Corr = AutoCorr(m_seq)
    Corr = zeros(1, numel(m_seq));
    for i = 1 : numel(m_seq)
        Corr(i) = sum(m_seq .* circshift(m_seq, i)) / numel(m_seq);
    end
end