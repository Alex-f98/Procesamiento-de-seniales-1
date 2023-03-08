
%% ACTIVIDAD 2: PASA ALTOS FLG.
% 
% 
% *1) Diseñar un filtro pasa altos que cumpla con las especificaciones indicadas 
% más abajo. Determine la ventana adecuada y el orden que cumplen con las especificaciones.* 
% 
% *2) Hallar la expresión de hd [n] e implementarla en Matlab. Graficar |H(𝜔)| 
% para 𝜔 ∈ [0,𝜋) y verificar si se satisfacen las especificaciones.* 
% 
% *3) Graficar h[n] y el diagrama de polos y ceros. ¿Qué tipo FLG resultó el 
% filtro implementado?* 
% 
% *4) Modificar el orden N en uno y repetir el punto 3). ¿Se siguen cumpliendo 
% las especificaciones?*
% 
% 
% 
% *0.992 < |H(w)| <1.008        0.7 < |w|<* $\pi$
% 
% *|H(w)| < 0.001             0 < |w| < 0.5*$\pi$

clc; clear; close all;
deltaP = 0.001; %0.05
dletaS = 0.008; %miro el mas chico, desde la de Hann.

DW = (0.7 - 0.5)*pi;
wc= (0.7 + 0.5)*pi/2;
%%
M= 12*pi/DW
N= M-1
n= 0:N;
%N = 28; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%iterar.
pb= @(nm, wc, N)( (wc/pi) .* sinc((wc/pi) .* (nm - N/2)) );
hd = pb(n, pi, N) - pb(n, pi-wc, N);

nfft = 8192;
%Blackman.
w = linspace(0, 2*pi, nfft);
w5 = window(@blackman, N+1)';

h = hd.*w5;
% FFT de ventanas
W5 = fft(h, nfft);

% Respuesta en el tiempo de las ventanas
figure
stem(w5), title('Black')%%% xlim([-29 N+30]), ylim([0 1.2])

figure
plot(w/pi, abs(W5));
xlim([0 1])
grid on;

%% 
% 
% 
% 
% 
%