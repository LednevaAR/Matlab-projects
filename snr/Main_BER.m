close all; clear; clc;
%% Init parametrs of model
Length_Bit_vector = 12e4;

Constellation = "16-QAM"; %_ 
Constellation2 = "BPSK"; %2
Constellation3 = "QPSK"; %3
Constellation4 = "8PSK"; %4
SNR = 30; % dB

%% Bit generator

Bit_Tx = randi([0,1], 1, Length_Bit_vector);

%% Mapping

IQ_TX = mapping(Bit_Tx, Constellation);

ForImagine = [0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 1 0 1 0 0 0 1 0 1 0 1 1 0 0 1 1 1 1 0 0 0 1 0 0 1 1 0 1 0 1 0 1 1 1 1 0 0 1 1 0 1 1 1 1 0 1 1 1 1];
ForImagine_Cons = mapping(ForImagine, Constellation);
scatterplot(ForImagine_Cons);   

dx = 0.05 + 0.05i; % displacement so the text does not overlay the data points
i = 0:15;
text(real(ForImagine_Cons(i+1)+dx), imag(ForImagine_Cons(i+1)+dx), dec2bin(i,4),'Color','white');
xlim([-1 1.4]);
ylim([-1 1.2]);
title(Constellation);

IQ_TX2 = mapping(Bit_Tx, Constellation2);
 
ForImagine = [0 1];
ForImagine_Cons = mapping(ForImagine, Constellation2);
scatterplot(ForImagine_Cons);   
i = 0:1;
text(real(ForImagine_Cons(i+1)+dx), imag(ForImagine_Cons(i+1)+dx), dec2bin(i,1),'Color','white');
xlim([-1 1.2]);
title(Constellation2);

IQ_TX3 = mapping(Bit_Tx, Constellation3);  

ForImagine = [0 0 0 1 1 0 1 1];
ForImagine_Cons = mapping(ForImagine, Constellation3);
scatterplot(ForImagine_Cons);   
i = 0:3;
text(real(ForImagine_Cons(i+1)+dx), imag(ForImagine_Cons(i+1)+dx), dec2bin(i,2),'Color','white');
xlim([-1 1]);
ylim([-1 1]);
title(Constellation3);

IQ_TX4 = mapping(Bit_Tx, Constellation4);

ForImagine = [0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1];
ForImagine_Cons = mapping(ForImagine, Constellation4);
scatterplot(ForImagine_Cons);   
i = 0:7;
text(real(ForImagine_Cons(i+1)+dx), imag(ForImagine_Cons(i+1)+dx), dec2bin(i,3),'Color','white');
xlim([-1 1.2]);
ylim([-1 1.2]);
title(Constellation4);


%% Channel
% Write your own function Eb_N0_convert(), which convert SNR to Eb/N0
Eb_N0 = Eb_N0_convert(SNR, Constellation);
% Use your own function of generating of AWGN from previous tasks
IQ_RX = Noise(IQ_TX, SNR);

%% Demapping
Bit_Rx = demapping(IQ_RX, Constellation);

%% Error check
% Write your own function Error_check() for calculation of BER
BER = Error_check(Bit_Tx, Bit_Rx);

%% Построим график BER(Eb/N0) от SNR для различных типов созвездий 

Len = 40;
counter = 1;
Ber = zeros(1, Len+10);
EbNo = zeros(1, Len+10);
Ber2 = zeros(1, Len+10);
EbNo2 = zeros(1, Len+10);
Ber3 = zeros(1, Len+10);
EbNo3 = zeros(1, Len+10);
Ber4 = zeros(1, Len+10);
EbNo4 = zeros(1, Len+10);
for i = -10:1:Len

    IQ_RX = Noise(IQ_TX, i);
    IQ_RX2 = Noise(IQ_TX2, i);
    IQ_RX3 = Noise(IQ_TX3, i);
    IQ_RX4 = Noise(IQ_TX4, i);

    Bit_Rx = demapping(IQ_RX, Constellation);
    Ber(counter) = Error_check(Bit_Tx, Bit_Rx);
    EbNo(counter) = Eb_N0_convert(i, Constellation);
    
    Bit_Rx2 = demapping(IQ_RX2, Constellation2);
    Ber2(counter) = Error_check(Bit_Tx, Bit_Rx2);
    EbNo2(counter) = Eb_N0_convert(i, Constellation2);
    
    Bit_Rx3 = demapping(IQ_RX3, Constellation3);
    Ber3(counter) = Error_check(Bit_Tx, Bit_Rx3);
    EbNo3(counter) = Eb_N0_convert(i, Constellation3);
    
    Bit_Rx4 = demapping(IQ_RX4, Constellation4);
    Ber4(counter) = Error_check(Bit_Tx, Bit_Rx4);
    EbNo4(counter) = Eb_N0_convert(i, Constellation4);

    counter = counter + 1;

end
%% 
p = figure;
semilogy(EbNo, Ber)
hold on 
semilogy(EbNo2, Ber2)
hold on 
semilogy(EbNo3, Ber3)
hold on 
semilogy(EbNo4, Ber4)

grid on
ylabel('BER')
xlabel('E_b/N_0 (dB)')

legend('16-QAM','BPSK','QPSK', '8PSK');
xlabel('Eb/N0 (dB)');
ylabel('BER');

title ("BER(Eb/N0)");
saveas(p, "BER(Eb_N0).fig");
m = figure;
semilogy(-10:1:Len, Ber)
hold on 
semilogy(-10:1:Len, Ber2)
hold on 
semilogy(-10:1:Len, Ber3)
hold on 
semilogy(-10:1:Len, Ber4)

legend('16-QAM','BPSK','QPSK', '8PSK');
xlabel(' SNR (dB)');
ylabel('BER');

title ("BER(SNR)");
saveas(m,"BER(SNR).fig")
%% Теоретический график
EbN = 5:1:22;
EbN0 = 10 .^ (EbN/10);
BER2 = (1/2)*erfc(sqrt(EbN0));
BER4 = (1/2)*erfc(sqrt(EbN0));
BER8 = (1/3)*erfc(sqrt(EbN0*3)*sin(pi/8));
BER16 = (1/2)*erfc(sqrt(EbN0*6/15));
f = figure;
semilogy(EbN, BER2)
hold on
semilogy(EbN, BER4)
hold on
semilogy(EbN, BER8)
hold on
semilogy(EbN, BER16)
hold off
grid on
title("Theoretical lines of BER(Eb/N0)")
legend('BPSK','QPSK', '8PSK','16-QAM')
ylabel('BER')
xlabel('Eb/N0(dB)')
saveas(f, "Theoretical.fig")
%% Additional task. Modulation error ration
MER_estimation = MER_my_func(IQ_RX, Constellation);

MerEr = zeros(1, 101);
MerEr2 = zeros(1, 101);
MerEr3 = zeros(1, 101);
MerEr4 = zeros(1, 101);
counter = 1;
for i = -50:1:50
    IQ_RX = Noise(IQ_TX, i);
    IQ_RX2 = Noise(IQ_TX2, i);
    IQ_RX3 = Noise(IQ_TX3, i);
    IQ_RX4 = Noise(IQ_TX4, i);
    
    MER_estimation = MER_my_func(IQ_RX, Constellation);
    MerEr(counter) = abs(MER_estimation - i);

    MER_estimation2 = MER_my_func(IQ_RX2, Constellation2);
    MerEr2(counter) = abs(MER_estimation2 -i );

    MER_estimation3 = MER_my_func(IQ_RX3, Constellation3);
    MerEr3(counter) = abs(MER_estimation3 - i);

    MER_estimation4 = MER_my_func(IQ_RX4, Constellation4);
    MerEr4(counter) = abs(MER_estimation4 -i);

    counter = counter + 1;
end

g = figure; 
plot ( -50:1:50, MerEr);
hold on 
plot ( -50:1:50, MerEr2);
plot ( -50:1:50, MerEr3);
plot ( -50:1:50, MerEr4);
hold off 
legend('16-QAM','BPSK','QPSK', '8PSK');
saveas(g, "Compare.fig")

%% Проверим, можно ли кривую BER(Eb/N0) с АБГШ  аппроксимировать
% теоретическая кривой.  
d = figure;
semilogy(EbNo, Ber, 'Color', 'green')
hold on 
semilogy(EbNo2, Ber2, 'Color', 'black')

semilogy(EbNo3, Ber3, 'Color', 'red')
 
semilogy(EbNo4, Ber4, 'Color', 'blue')

semilogy(EbN, BER2,  'Color', 'black')

semilogy(EbN, BER4, 'Color', 'red')

semilogy(EbN, BER8, 'Color', 'blue')

semilogy(EbN, BER16,'Color', 'green')
hold off

grid on
ylabel('BER')
xlabel('E_b/N_0 (dB)')

legend('16-QAM','BPSK','QPSK', '8PSK');
xlabel('Eb/N0 (dB)');
ylabel('BER');
xlim ([5 10]);
ylim([ 8e-6 0.5]);
title ("Сравнение результатов");
saveas(d, "CompareTeoryAndPractice.fig")

% Compare the SNR and MER_estimation from -50dB to +50dB for BPSK, QPSK,
% 8PSK and 16QAM constellation.
% Plot the function of error between SNR and MER for each constellation 
% Discribe the results. Make an conclusion about MER.
% You can use the cycle for collecting of data
% Save figure

%% Experimental BER(SNR) and BER(Eb/N0)
% Collect enough data to plot BER(SNR) and BER(Eb/N0) for each
% constellation.
% Compare the constalation. Describe the results
% You can use the cycle for collecting of data
% Save figure

%% Чем больше точек в созвездии, тем вероятнее ошибка при одинковом SNR

%% Theoretical lines of BER(Eb/N0)
% Read about function erfc(x) or similar
% Configure the function and get the theoretical lines of BER(Eb/N0)
% Compare the experimental BER(Eb/N0) and theoretical for BPSK, QPSK, 8PSK
% and 16QAM constellation
% Save figure