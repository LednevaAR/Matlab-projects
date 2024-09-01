clear; clc; close all;
%% 
Header = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%Header = [1 1 1 1 1 1 1 1];
%%
Stream = importdata('Matlab_L3_3.mat');
Frame_Corr = zeros(1, numel(Stream) - numel(Header));
Frame = Frame_Corr;
Pred = 0;
for i = 1 : (numel(Stream) - numel(Header))
    Frame_Corr(i) = sum(Header .* Stream((i + 1) : (i + numel(Header))))/numel(Header);
    if (Frame_Corr(i) == 1)
        Diff = i - Pred;
        if ((Diff >= 884) || (Pred == 0))
            Pred = i;
            Frame(i) = Frame_Corr(i);
        end
    end
end
%%
Start_Of_Frame_Position = find(Frame_Corr == 1, 1);
%%
Number_Of_Frame = numel(find(Frame > 0));
F = find(Frame > 0);
D = zeros(1, numel(F) - 1);
for i = 1:(numel(find(Frame > 0)) - 1)
    D(i) = F(i + 1) - F(i);
end
L = round(mean(D));
disp(find(Frame > 0));
%%
f = figure();
xlabel('Bit number');
title('Correlation graph');
plot(Frame_Corr);
saveas(f, 'Frame_Corr.fig');
%%
Frame_search = [Start_Of_Frame_Position, Number_Of_Frame];
save('Frame_search.mat', 'Frame_search');
%%
Length_of_Data = (L - numel(Header)) * Number_Of_Frame;
fprintf('Length_of_Data = %d', Length_of_Data);