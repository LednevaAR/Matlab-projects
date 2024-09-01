% Clear workspace. Clear command window. Close all figure
clear; clc; close all;
%% task
% 0) Create a function
% 1) Create the random matrix
% 2) Analyse the code. Insert the calculation of runtime of code
% 3) rewrite the code in more optimised form in matlab
% 4) Provide the evidence that results matrix and legacy matrix is the same
% 5) calculate the runtime of new code. Compare it with legacy code. Make
% an conclusion about code. Which one is the more optimised? Which code do
% you suggest to use in matlab? And why?
%% Config the system

% Fixed random generator
rng(11);
% TO-DO 1%
% create Input_Matrix matrix 20-to-750 size and
% with normal distributed numbers from -33 to 100
%
Input_Matrix = Matrix_generator();
%disp(Input_Matrix);
Legacy_Matrix = Input_Matrix;
Ethalon_Matrix = Input_Matrix;
%% Run legacy code
% TO-DO 2
% Measure the runtime of current function
T = zeros(1, 100);
for i = 1:1000
    tStart = tic;
    Legacy_Output_Matrix = Legacy_Instruction(Input_Matrix);
    T(i) = toc(tStart);
end
Time_legacy_code = mean(T);
% Save the runtime in variable
% Time_legacy_code = TIME;

%% Run optimised code
% TO-DO 3
% Measure the runtime of your function
% Create function New_Instruction()
% Rewrite and optimised function Legacy_Instruction()
% Use matrix operation, logical indexing
% Try not to use the cycle
for i = 1:1000
    tStart = tic;
    Optimised_Output_Matrix = New_Instruction(Input_Matrix);
    T(i) = toc(tStart);
end
Time_Optimised_code = mean(T);
% Save the runtime in variable
% Time_Optimised_code = TIME;
    
%% Checking the work of student
% TO-DO 4
% Compare the matrix and elapsed time for instruction
% Matrix must be equal each other, but the runtime sill be different

% Runtime comparison
% Comparison of matrix
% Matrix size and value
if (size(Legacy_Output_Matrix) == size(Optimised_Output_Matrix))
    if (size(find((Legacy_Output_Matrix - Optimised_Output_Matrix) ~= 0), 1) == 0)
        fprintf('Time_legacy_code = %f\n', Time_legacy_code);
        fprintf('Time_Optimised_code = %f\n', Time_Optimised_code);
        disp("Legacy_Output_Matrix and Optimised_Output_Matrix are same.")
        if Time_Optimised_code < Time_legacy_code
            disp('Optimised function is faster than legacy.');
        else
            disp('Legacy function is faster than optimised.');
        end
        fprintf("Size of matrix = %dx%d\n", size(Legacy_Output_Matrix, 1), size(Optimised_Output_Matrix, 2));
        fprintf('Conclusion: Optimised function is %f times faster than legacy.\nThis happens because using indexing is faster and more productive than using loops.\n', Time_legacy_code/Time_Optimised_code);
    else
        disp("Legacy_Output_Matrix and Optimised_Output_Matrix are different.")
    end
else
    disp("Size of Legacy_Output_Matrix and size of Optimised_Output_Matrix are different.")
end

        


%% Function discribing

function Output_Matrix = Legacy_Instruction(Matrix)
   
    for itter_rows = 1 : size(Matrix,1)
        for itter_column = 1 : size(Matrix,2)
            if mod(itter_rows,3) == 0
                Matrix(itter_rows,itter_column) = -Matrix(itter_rows,itter_column);
            end
        end
    end

   for itter_rows = 1 : size(Matrix,1)
        for itter_column = 1 : size(Matrix,2)
            if Matrix(itter_rows,itter_column) >= 0
                Matrix(itter_rows,itter_column) = Matrix(itter_rows,itter_column) + 55.2;
            end
        end
   end

   for itter_rows = 1 : size(Matrix,1)
        for itter_column = 1 : size(Matrix,2)
            if Matrix(itter_rows,itter_column) > 35
                Matrix(itter_rows,itter_column) = Matrix(itter_rows,itter_column)^2;
            end
        end
    end

    Output_Matrix = Matrix;
end

function New_Matrix = Matrix_generator()
    New_Matrix = rand(20, 750) * 133 - 33;
end

function Output_Matrix = New_Instruction(Matrix)
    Matrix(3:3:end, :) = -Matrix(3:3:end, :);
    Matrix = Matrix + (Matrix >= 0)*55.2;
    Matrix = Matrix .* (Matrix <= 35) + Matrix .* Matrix .* (Matrix > 35);
    Output_Matrix = Matrix;
end