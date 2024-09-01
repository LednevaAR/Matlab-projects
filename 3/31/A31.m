clear; clc; close all;
%%
[audio, Fs] = audioread("file4.wav");
Spec = fft(audio);
f = (0:(Fs/length(audio)):length(audio)/2);
Cell = cell(numel(f) - 1, 2);
for i = 0:(numel(f) - 2)
    Cell{i + 1, 1} = f(i + 1);
    Cell{i + 1, 2} = abs(Spec(i + numel(f)));
end
figure;
n = length(audio);
fr = (-n/2:n/2-1)*(Fs/n);
plot(fr, abs(Spec))
xlabel('Frequency');
ylabel('Amp');
title('Spec');
Max = sortrows(Cell, 2, 'descend');
fprintf('Harmonic frequencies = %d Hz, %d Hz and %d Hz', Max{1, 1}, Max{2, 1}, Max{3, 1});
%%
FilteredSpec = Filter(Spec, Max(1:3, 1));
figure;
n = length(audio);
fr = (-n/2:n/2-1)*(Fs/n);
plot(fr, abs(FilteredSpec))
title('FilteredSpec')
xlabel('Frequency');
ylabel('Amp');
new_audio = ifft(FilteredSpec);
new_audio = new_audio/max(abs(new_audio));
audiowrite('result4.wav', new_audio, Fs)