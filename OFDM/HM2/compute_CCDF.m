function [CCDF, PAPR_sorted] = compute_CCDF(PAPR_slid_wind)
    PAPR_sorted = sort(PAPR_slid_wind);
    CCDF = zeros(1, numel(PAPR_sorted));
    for i = 1 : length(PAPR_sorted)
        CCDF(i) = sum(PAPR_sorted > PAPR_sorted(i)) / numel(PAPR_slid_wind);
    end
end

