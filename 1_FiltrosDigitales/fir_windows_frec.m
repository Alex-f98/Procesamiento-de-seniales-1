clear all
close all

%% Ventanas para diseÃ±o de filtros FIR

% Orden del filtro
N = 200;

% Vector de frecuencias
nfft = 8192;
w = linspace(0,2*pi,nfft);

% Ventanas con ripple fijo
w1 = window(@rectwin, N+1)';
w2 = window(@bartlett, N+1)';
w3 = window(@hamming, N+1)';
w4 = window(@hann, N+1)';
w5 = window(@blackman, N+1)';

% FFT de ventanas
W1 = fft(w1,nfft);
W2 = fft(w2,nfft);
W3 = fft(w3,nfft);
W4 = fft(w4,nfft);
W5 = fft(w5,nfft);

% Normalización para comparar
W1 = W1/W1(1);
W2 = W2/W2(1);
W3 = W3/W3(1);
W4 = W4/W4(1);
W5 = W5/W5(1);


%% Graficos

% Respuesta en frecuencia de las ventanas
figure
subplot(5,10,[1 2 3])
plot(w/pi,20*log10(abs(W1)),'Color',[0, 113, 188]/256,'LineWidth',2), ylim([-100 0])
ylabel(['Magnitud [dB]';'  (Rectwin)  ']), xlabel('\omega / \pi'), xlim([0 0.05]), grid on

subplot(5,10,[11 12 13])
plot(w/pi,20*log10(abs(W2)),'color',[236, 176, 31]/256, 'LineWidth',2), ylim([-100 0])
ylabel(['Magnitud [dB]';'  (Bartlett) ']), xlabel('\omega / \pi'), xlim([0 0.05]), grid on

subplot(5,10,[31 32 33])
plot(w/pi,20*log10(abs(W3)),'color',[216, 82, 24]/256,'LineWidth',2), ylim([-100 0])
ylabel(['Magnitud [dB]';'  (Hamming)  ']), xlabel('\omega / \pi'), xlim([0 0.05]), grid on

subplot(5,10,[21 22 23]), 
plot(w/pi,20*log10(abs(W4)),'color',[125, 46, 141]/256,'LineWidth',2), ylim([-100 0])
ylabel(['Magnitud [dB]';'    (Hann)   ']), xlabel('\omega / \pi'), xlim([0 0.05]), grid on

subplot(5,10,[41 42 43])
plot(w/pi,20*log10(abs(W5)),'color',[118, 171, 47]/256,'LineWidth',2), ylim([-100 0])
ylabel(['Magnitud [dB]';'  (Blackman) ']), xlabel('\omega / \pi'), xlim([0 0.05]), grid on

subplot(5,10,[[5 6 7 8 9 10]  10+[5 6 7 8 9 10]  20+[5 6 7 8 9 10] 30+[5 6 7 8 9 10] 40+[5 6 7 8 9 10]])
plot(w/pi,db([abs(W1)' abs(W3)'  abs(W2)' abs(W4)' abs(W5)']),'LineWidth',2)
legend('Rectwin','Hamming','Bartlett','Hann', 'Blackman')
ylabel('Magnitud [dB]')
xlabel('\omega / \pi (Frecuencia normalizada)')
ylim([-150 20])
xlim([0 1])
grid on
title(['N=',num2str(N)])
set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % full screen. 

