clear; clc; close all;
%% 1
[alphabet, emp_prob] = alphabet_probabilities('aaagggghhhhh');
disp('Алфавит')
disp(alphabet)
disp('Эмпирические вероятности')
disp(emp_prob)
%% 2
H = entropy('aaabbbccc');
disp('Безусловная энтропия равна')
disp(H)
%% 3
load strings
H = cond_val_entropy(X, Y, 'a');
if (H ~= -1)
    sprintf('H(X|Y = y) = %g', H)
else
    sprintf('В Y нет символа %c', y)
end
%% 4
load strings
H1 = cond_entropy(X, Y);
H2 = joint_entropy(X, Y);
H3 = cond_entropy(X, Y) + entropy(Y);
H4 = joint_entropy1(X, Y);
sprintf('H(X|Y) = %g', H1)
sprintf('H(X,Y) = %g', H2)
disp(H3)
disp(H4)
%% 5
load strings
disp('-----------------------------------------------------------------------------')
HX = entropy(X);
sprintf('H(X) = %g', HX)
HY = entropy(Y);
sprintf('H(Y) = %g', HY)
[N, ~] = alphabet_probabilities(X);
RX = 1 - (entropy(X) / log2(length(N)));
sprintf('R(X) = %g', RX)
[M, ~] = alphabet_probabilities(Y);
RY = 1 - (entropy(Y) / log2(length(M)));
sprintf('R(Y) = %g', RY)
H_YX = cond_entropy(Y, X);
sprintf('H(Y|X) = %g', H_YX)
H_XY = cond_entropy(X, Y);
sprintf('H(X|Y) = %g', H_XY)
HXY = joint_entropy(X, Y);
sprintf('H(X,Y) = %g', HXY)
IXY = HX - H_XY;
sprintf('I(X;Y) = %g', IXY)
IYX = HY - H_YX;
sprintf('I(Y;X) = %g', IYX)
HXY1 = HX + HY - IXY;
sprintf('H(XY) = %g', HXY1)
if ((H_XY < HX) && (HX < log2(length(N))))
    disp('H(X|Y) <= H(X) <= log2(L)')
end
if ((H_YX < HY) && (HY < log2(length(M))))
    disp('H(Y|X) <= H(Y) <= log2(L)')
end
if (round(IXY - IYX) == 0)
    disp('I(X;Y) = I(Y;X)')
end
if (HXY1 - (HX + H_YX) == 0)
    disp('H(X,Y) = H(X) + H(Y|X)')
end
if (HXY1 - (HY + H_XY) == 0)
    disp('H(X,Y) = H(Y) + H(X|Y)')
end
%% 6 1
load strings
disp('-----------------------------------------------------------------------------')
L = 2;
[H_C, H_L] = f(L, X);
sprintf('H_C(L) = %g', H_C)
sprintf('H_L = %g', H_L)
[A, ~] = alphabet_probabilities(X);
sprintf('Энтропия нулевого порядка равна %g', log2(length(A)))
sprintf('Энтропия символа равна %g', entropy(X))
n = length(X);
num = 2 : n;
HC = zeros(1, n - 1);
HL = zeros(1, n - 1);
for i = 1 : n - 1
    [HC(i), HL(i)] = f(i + 1, X);
end
figure
plot(num, HC, "-o", num, HL, "-o", 'MarkerSize', 2)
legend('H_C', 'H_L')