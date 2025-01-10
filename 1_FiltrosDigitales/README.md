# Procesamiento de señales I (86.51) 
# Trabajo Practico I: Diseño de filtros diigitales

## 1. Introducción

Objetivo Se requiere procesar una señal de audio transmitida entre dos equipos a través de un determinado medio. En el equipo receptor, la señal se recibe contaminada por varias interferencias de banda angosta, como se observa en la Figura 1.

![Filtrado de una señal acustica contaminada con interferencias](./img/imagen1.png)


Para mitigar este problema, afectando lo menos posible a la señal original, se propone el diseño e implementación de un filtro FIR digital de Fase Lineal Generalizada. El diseño deberá contemplar múltiples eliminabandas selectivos en frecuencia.
### 1.2.
 Requerimientos La selectividad deberá ser definida de manera tal que el ancho de banda utilizado para suprimir cada componente de interferencia no exceda los 300 Hz. Se requiere que la máxima variación en las bandas de paso y suprimida no superen δp = 0, 08 y δs = 0, 016 respectivamente.
## 2. Desarrollo
*Problema 1* 
(a) 

Utilice la función audioread() de Matlab para cargar la señal de interferencias. Grafíque el espectro de dicha señal y verifique las frecuencias de las tres componentes principales

```matlab
clear; close all; clc
addpath("C:.\FiltrosDigitales\TP1_FiltrosDigitales\CANCIONES");

%Señales
[S1, Fs_1] = audioread("Pista_01.wav");
[S2, Fs_2] = audioread("Pista_02.wav");
[S3, Fs_3] = audioread("Pista_03.wav");
[S4, Fs_4] = audioread("Pista_04.wav");
[S5, Fs_5] = audioread("Pista_05.wav");
%interferencias
[V1, Fs_v1] = audioread("interferencias.wav");
%[V2, Fs_v2] = audioread("interferencias2.wav");
f1 = 1400; % Hz
f2 = 2735; % Hz
f3 = 3772; % Hz
```
Para definir los puntos se toma la longitud de puntos del ruido, pero estas son demasiadas y MATLAB (al menos en LiveScript) tiende a recortar los tamaños para reducir tiempo, así que se lo achica pero aun así es indistinto.
Se procede a graficar el espectro del ruido en escala logarítmica para confirmas las frecuencias dadas

```matlab
NFFT= length(V1)/10;
% Ejes de frecuencia.
F = ( (0 : 1/NFFT : 1-1/NFFT) * Fs_v1).';
w = linspace(0, 2*pi, NFFT);

figure(1)
    FFTN = fft(V1, NFFT);  %Transformada de fourier
    plot( F(1:NFFT/2), 20*log10( abs( FFTN(1:NFFT/2) ) ), 'r');
    xlabel("Frecuencia [Hz]")
    ylabel("Magnitud [dB]")
    xlim([0 5000])
    title("Trasformada discreta de fourier del ruido")
    %f1 = 1400 Hz, f2 = 2735Hz y f3 = 3772 Hz
    
    hold on
    plot([1400, 1400], [-80,80],'--k')
    plot([2735, 2735], [-80,80],'--k')
    plot([3772, 3772], [-80,80],'--k')
    grid on
```

