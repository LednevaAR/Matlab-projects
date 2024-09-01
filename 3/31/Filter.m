function Filtered = Filter(Spec, Max)
    F = zeros(numel(Spec), 1);
    for i = 1:3
        a = Max{i, 1} + round(numel(Spec)/2) + 1;
        F(a) = 1;
        a = round(numel(Spec)/2) - Max{i, 1} + 1;
        F(a) = 1;
    end
    Filtered = Spec .* F;
end
