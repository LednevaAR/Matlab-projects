clear; clc; close all;
%%
Register = [1 0 0 1 0 1 0 1 0 0 0 0];
PN_seq = RSLOS(Register);
%PN_seq = randi([0, 1], 1, 1e3);
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
%%
PeriodPN = zeros(1, numel(Register));
for i = 1 : numel(Register)
    PN_seq = RSLOS(Register(1:i));
    AutoCorrPN = AutoCorr(PN_seq);
    [M, PeriodPN(i)] = max(AutoCorrPN(1:(numel(AutoCorrPN) - 1)));
end
f = figure();
plot(PeriodPN);
xlabel('Register length');
ylabel('Period of PN');
title('PeriodPN function');
grid on;
saveas(f, 'PN_Period.fig')
%%
saveas(f, 'PN_Period.png')
