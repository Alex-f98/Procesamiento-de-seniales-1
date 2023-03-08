clear all
close all

%% Ventanas para diseño de filtros FIR

% Orden del filtro
N = 200;

% Ventanas con ripple fijo
w1 = window(@rectwin, N+1)';
w2 = window(@bartlett, N+1)';
w3 = window(@hamming, N+1)';
w4 = window(@hann, N+1)';
w5 = window(@blackman, N+1)';


%% Graficos

% Respuesta en el tiempo de las ventanas
figure
subplot(3,2,[1 2])
stem(w1), title('Rectwin'), xlim([-29 N+30]), ylim([0 1.2])
subplot(3,2,3)
stem(w2), xlim([1 N+1]), title('Bartlett')
subplot(3,2,4)
stem(w3), xlim([1 N+1]), title('Hamming')
subplot(3,2,5)
stem(w4), xlim([1 N+1]), title('Hann')
subplot(3,2,6)
stem(w5), xlim([1 N+1]), title('Blackman')
set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % full screen. 


