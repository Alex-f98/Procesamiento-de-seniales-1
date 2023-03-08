% Defina la señal x(n) captada por MIC-1 utilizando como señal de interés s(n) cualquiera
% de los archivos de audio de muestra disponibles para el TP1 y como ruido ambiente v(n) a
% una secuencia de ruido blanco gaussiando de media nula y varianza σ_v^2
% tal que se cumpla una relación señal a ruido sea SNR = 10 dB. Luego defina la señal captada por MIC-2 como un
% proceso de media móvil generado a partir de v(n): u(n) = 0, 8 v(n)+ 0, 2 v(n−1)−0, 1 v(n−2).
clc; clear; close all
%MIC-1: x(n) = s(n) + v(n) / s(n) del TP1(CANCIONES), v(n) RBG(mu=0, varR=sigmav2)
% tal que SNR = 10dB
[s_n, fs] = audioread('CANCIONES/Pista_01.wav');     Ts = 1/fs;
sn = s_n(1:20e4);
L = length(sn);                                     n = 0:L-1;
%varS = var(sn);
SNR_dB = 10;                                        %SNR = 10log(varS/varR)

[xn, un, vn, varR] = fun(sn, SNR_dB);

%% a) Filtrado: solo se dispone de la señal observada x(n)(MIC-1) y el proceso u(n)(MIC-2)
mu_LMS = 0.5;
M = 2;
w_inicial = [3; 4];

%[Weigths, errors, x_est] = LMS(un, xn, M, mu_LMS, w_inicial);
[Weigths, errors, e_V, x_est] = LMS2(un, xn, vn, M, mu_LMS, w_inicial); 
%Una sola realizacion:
Vn_ = e_V.^2; % abs(vn(1:end-1)' - x_est).^2;

figure()
    plot(Vn_, '-b', 'LineWidth', 2)
    grid on;
    title("Diferencia al cuadrado entre el ruido y la salida del filtro")
    ylabel("$\hat{V}(n) = | v(n) -\hat{d}(n) |^2$", 'Interpreter', 'latex')
    xlabel("Nro de iteracion")
%%
%Grafico los coeficientes    
figure()
    hold on
    plot(Weigths(1,:), 'r', 'LineWidth', 2)
    plot(Weigths(2,:), 'b', 'LineWidth', 2)
    
    title("Coeficientes estimados")
    ylabel("w(n)", 'Interpreter', 'latex')
    xlabel("Nro de iteracion", 'Interpreter', 'latex')
    grid on

%-% Grafico wn_ en funcion de las iteraciones en el plano (w1 w2), curva de nivel J(w)
CurvaLvL_J(un, xn, Weigths, M, L)
title("Trayectoria de los pesos, con w(0)=[3,4]^T")

%% b)Reproducción

%1. Pista de audio contaminada, x(n).
    disp("1. Pista de audio contaminada, x(n).");
    sound(xn, fs); pause(length(xn)/fs); clear sound
    %Ruido bastante apreciable
    
%2. Señal de error, e(n).
    disp("2. Señal de error, e(n).");
    sound(errors, fs); pause(length(errors)/fs); clear sound
    %Ruido apenas apreciable, pero se hace presente
    
%3. Pista de audio original, s(n).
    disp("3. Pista de audio original, s(n).");
    sound(sn, fs); pause(length(sn)/fs); clear sound
    %Señal original, sin ruido, el constraste con el error es apreciable
    %prestando atencion.
     
   
    
    
    
%%




% function [xn, un, vn, varR] = fun(sn, SNR_dB)
%     L = length(sn);
%     varS = var(sn);
%     varR = varS / 10^(SNR_dB/10);  
%     vn = sqrt(varR)*randn(L, 1);
% 
%     xn = sn + vn;
% 
%     %MIC-2: v(n): u(n) = 0.8 v(n)+ 0.2 v(n−1)−0.1 v(n−2).
%     h = [0.8 0.2 -0.1];
%     un = filter(h, 1, vn);
% 
% end


