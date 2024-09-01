clear; clc; close all;

%% Analisys of text from file
% to-do
% 1) Reading the file as a text of array of chars
% 2) Create array of cells which consist of three columns
% "char"->"amount of meetings"->"probobilities"
% 3) Save the chars and probobalities to file *.mat and *.xls as the cell
% variables. Name the files shouold be:
% Data_Analisys.mat
% Data_Analisys.xls
% Data_Analisys.png
% 4) Plot the distribution of probability of symbols in text. 
% Be careful to the labels on the axis.
% Recommendation use xticks(), xticklabels().
% 5) Save the plot as figure and PNG image with resolution at least 400 px. The name
% of files should be: Data_Analisys.png
%% Reading the file
% TO-DO 1
% Read the file from *.txt as a char stream
%fid = fopen("Textvar3.txt", "r");
fid = fopen("Text1.txt", "r");
symbols = lower(fread(fid, '*char'));

%% Analysis
% TO-DO 2 
% Use ony char from file
% Use lowercase string
% Try to use the "Cell" as a data containers;
% Name the varible Data_Analisys
% The cell should consist of 3 columns:
% "Symbol" | "Amount of meeting" | "Probolitie"

% You can use only 1 cycle for this task
% Avoid the memmory allocation in cycle
Sort = unique(symbols);
Data_Analisys = cell(numel(Sort), 3);
for i = 1:numel(Sort)
    Data_Analisys{i, 1} = Sort(i);
    Data_Analisys{i, 2} = numel(find(symbols == Sort(i)));
    Data_Analisys{i, 3} = numel(find(symbols == Sort(i)))/numel(symbols);
end
DA = sortrows(Data_Analisys, 3, 'descend');
%% Plot Data
% TO-DO 3
% Illustrate the results from Analysis block
% There should be lable of axis, title, grid
f = figure();
f.WindowState = 'maximized';
X = linspace(1, size(Data_Analisys, 1), size(Data_Analisys, 1));
Y = cell2mat(Data_Analisys(:, 3));
bar(X, Y, 'g');
xlabel('Symbol');
ylabel('Probability');
title('Probability of each symbol in text');
grid("on");
ylim([0, 0.14]);
xticks(X);
Data = string(Data_Analisys(:, 1));
N = find(Data == newline);
R = find(Data == char(13));
Data(N(1)) = '\n';
Data(R(1)) = '\r';
xticklabels(Data);
f1 = figure;
X = linspace(1, size(DA, 1), size(DA, 1));
Y = cell2mat(DA(:, 3));
bar(X, Y, 'g');
xlabel('Symbol');
ylabel('Probability');
title('Probability of each symbol in text');
xticks(X);
Data = string(DA(:, 1));
N = find(Data == newline);
R = find(Data == char(13));
Data(N(1)) = '\n';
Data(R(1)) = '\r';
xticklabels(Data);

%% Save the file
% TO-DO 4
% Save the figure as Data_Analisys.fig
saveas(f, 'Data_Analisys.fig')
% Save the figure as image Data_Analisys.png
saveas(f,'Data_Analisys.png')
saveas(f1,'DA.png')
% Save the data as MAT-file Data_Analisys.mat
save('Data_Analisys.mat', 'Data_Analisys')
% Save the data as Excel table Data_Analisys.xls
writecell(Data_Analisys, 'Data_Analisys.xls')

fclose('all');