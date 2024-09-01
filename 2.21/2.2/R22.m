clear; clc; close all;
%%
Register = [1 0 0 1 0 1 0 1 0 0 0 0];
PN_seq = RSLOS(Register);
%%
AutoCorrPN = AutoCorr(PN_seq);
BN = linspace(1, numel(AutoCorrPN), numel(AutoCorrPN));
f = figure();
plot(AutoCorrPN);
xlabel('Bit number');
ylabel('AutoCorrelation function');
title('AutoCorrelation function of PN-seq');
saveas(f, 'ACF_Srambler.fig')
%%
[M, Period] = max(AutoCorrPN(1:(numel(AutoCorrPN) - 1)));
fprintf('Generator_Period = %d\n', Period);
PN_Period = Period;
fprintf('PN_Period = %d\n', PN_Period);
%% проверка того, что период последовательности действительно равен Generator_Period
for i = 1 : numel(PN_seq)
    if isequal(PN_seq(1 : (i * floor(numel(PN_seq)/i))), repmat(PN_seq(1 : i), 1, floor(numel(PN_seq)/i))) && isequal(PN_seq((i * floor(numel(PN_seq)/i) + 1):end), PN_seq(1 : mod(numel(PN_seq), i)))
        PN_Period = i;
        break;
    end
end