function Corr = AutoCorr(PN)
    Corr = zeros(1, numel(PN));
    for i = 1 : numel(PN)
        Corr(i) = sum(PN .* circshift(PN, i))/numel(PN);
    end
end