%Bancos de filtros

%cheqear a la salida la Y constante y ver cuanto vale el retardo de todo
%completo.

%% e)BANCO DE FILTRO QMF.
% e) Implemente en Matlab el banco de filtros para el caso QMF. Suponiendo que se aplica un
% impulso como entrada x(n) = 𝛿(n), para c0 = c1= 1/√2 y n0 = n1 = 0, obtenga la salida y(n),
% grafiquela en el tiempo y en frecuencia. ¿Se alcanza la condición de PR?. ¿Cuánto es el
% retardo del sistema completo?

%h0(n) = c0.delta(n - 2n0) + c1.delta(n - 2n1 - 1)
%h1(n) = c0.delta(n - 2n0) - c1.delta(n - 2n1 - 1)
%f0(n) = c0.delta(n - 2n0) + c1.delta(n - 2n1 - 1)
%f1(n) = -c0.delta(n - 2n0) + c1.delta(n - 2n1 - 1)

%x = delta(n);
clc; clear; close all;
nfft = 10e3;
w = linspace(0, 2*pi, nfft);

c0 = 1/sqrt(2);
c1 = c0;
n0 = 0; 
n1 = n0;

%Entrada: delta(n);

x = [1, zeros(1,10)]; %[1 0 0 ... 0]

h00= c0*[zeros(1, 2*n0) 1 zeros(1, 2*n1+1)];
h01= c1*[zeros(1, 2*n1+1) 1 zeros(1, 2*n0)];

h0 = h00 + h01;
n = 0:length(h0)-1;

f0 = h0;
h1 = (-1).^n.*h0; %exp(i*pi*n) dexplazamiento en frecuencia == -h01
M = 2;

f1 = -h1;

%x->|H0| ->M---->L->|F0|--y0
%                          (+)-->yk                     
%x->|H1| ->M---->L->|F1|--y1
y0 = filter(f0, 1, upsample( downsample( filter(h0, 1, x), M ), M) );
y1 = filter(f1, 1, upsample( downsample( filter(h1, 1, x), M ), M) );

yk = [y0; y1];

y_n = sum(yk);
Y_n = fft(y_n, nfft);
%Se alcanza reconstruccion perfecta, el retardo es de 1
figure()
    hold on
    stem(0:length(x)-1, x, "LineWidth", 2)
    stem(0:length(y_n)-1, y_n, "LineWidth", 2)
    title("Salida del QMF (x_{input} = [1 0 0..0 0])")
    xlabel("n")
    ylabel("x(n),y(n)")
    legend("x(n)_{input}", "y(n)_{output}")
    grid on
    
%f)
% Grafique superpuestas las respuestas en frecuencia de |H0(𝜔)| y |H1(𝜔)|. Observe las
% transiciones y bandas de supresión de ambos filtros. También grafique la salida para cada
% rama por separado (Y0(w) e Y1(w)).
% Ayuda: para obtener la salida de cada rama puede anular la entrada sobre la otra rama.

figure()
    hold on;
    plot(w/pi, abs(Y_n));
    plot(w/pi, abs(fft(h0,nfft)))
    plot(w/pi, abs(fft(h1,nfft)))
    plot(w/pi, abs(fft(y0,nfft)))
    plot(w/pi, abs(fft(y1,nfft)))
    legend("Y_n", "H0", "H1", "Y0", "Y1")
    grid on;
    xlabel("w/ \pi")
    
%% 2)Bancos de filtro: Arbol.
clear; close all; clc;
% Suponga un banco de filtros de cuatro sub-bandas con una estructura de árbol 
% cómo el de la Figura. Asuma que posee filtros QMF: H0(z) = (1+z-1)/√2
%a) 
% Implemente el banco de filtros completo, banco de análisis y síntesis, 
% conectados directamente en cada sub-banda. Gráfique la salida y(n) para 
% una entrada impulsiva y su respuesta en frecuencia |Y(𝜔)|. 
% Se cumple la reconstrucción perfecta?
 M=2;
 L=2;
% H0 = (1 - z^-1)/sqrt(2);
h0 = [1 -1]./sqrt(2);       
x = [1 zeros(1, 10)];

%up
[x0, x1] = QMF_up(x, h0, M);            

[x00, x01] = QMF_up(x0, h0, M);         
[x10, x11] = QMF_up(x1, h0, M);         

%down
v00 = x00; v01 = x01; v10 = x10; v11 = x11;

v0 = QMF_down(v00, v01, h0, L);         
v1 = QMF_down(v10, v11, h0, L);         

y = QMF_down(v0, v1, h0, L);            

figure()
    hold on
    stem(0:length(x)-1, x, "LineWidth", 2)
    stem(0:length(y)-1, y, "LineWidth", 2)
    title("Salida del QMF-arbol (x_{input} = [1 0 0..0 0])")
    xlabel("n")
    ylabel("x(n),y(n)")
    legend("x(n)_{input}", "y(n)_{output}")
    grid on
    
 % retardo de un QMF es 1 (Nqmf/2) / Nqmf ={ [(N/2)*2]+ N }/2 = N
 %entonces el de este QMF-arbol es:
 %  N_QMF-arbol ={ [(N/2 + N)*2]+N }/2 = { 4N }/2
 % si N = 1;
 
 
%% 2) Banco de filtros: octava
clc; clear; close all;
%Suponga un banco de filtros de cuatro sub-bandas en octavas, ver Figura. Asuma que cada etapa posee filtros QMF: H0(z) = (1+z-1)/√2.
%a)
% Calcule los retardos D1 y D2 para compensar el delay de diferencia que se
% produce por los retardos introducidos en los filtros de las distintas sub-bandas.
%b) 
% Implemente el banco de filtros completo y grafique la respuesta en frecuencia de cada una de
% las cuatro ramas y la respuesta completa (es decir con las cuatro sub-bandas conectadas).

M=2;
L=2;
% H0 = (1 - z^-1)/sqrt(2);
h0 = [1 -1]./sqrt(2);       
x = [1 zeros(1, 10)];
n_D2 = 9;
n_D1 = 3;

%up
[x0, x1] = QMF_up(x, h0, M);
v1 = delay(x1, n_D2);

[x00, x01] = QMF_up(x0, h0, M);         
v01 = delay(x01, n_D1);

[x000, x001] = QMF_up(x00, h0, M); 
v000 = x000; v001 = x001;
%down


v00 = QMF_down(v000, v001, h0, L); 

v0 = QMF_down(v00, v01, h0, L);         

y = QMF_down(v0, v1, h0, L);            

figure()
    hold on
    stem(0:length(x)-1, x, "LineWidth", 2)
    stem(0:length(y)-1, y, "LineWidth", 2)
    title("Salida del QMF-arbol (x_{input} = [1 0 0..0 0])")
    xlabel("n")
    ylabel("x(n),y(n)")
    legend("x(n)_{input}", "y(n)_{output}")
    grid on
    

 %PREGUNTAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 %PREGNTAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
%%
function [x0, x1] = QMF_up(x, h0, M)
    %c0 = 1/sqrt(2);
%     c1 = c0;
% 
%     %n0 = 0; 
%     n1 = n0;
% 
%     h00= c0*[zeros(1, 2*n0) 1 zeros(1, 2*n1+1)];
%     h01= c1*[zeros(1, 2*n1+1) 1 zeros(1, 2*n0)];
% 
%     h0 = h00 + h01;
%     
     n = 0:length(h0)-1;
     h1 = (-1).^n.*h0; %exp(i*pi*n) dexplazamiento en frecuencia 
    
%     
    x0 = downsample(conv(h0, x), M);
    x1 = downsample(conv(h1, x), M);
end

function v0 = QMF_down(v00, v01, h0, L)
    
    n = 0:length(h0)-1;
    h1 = (-1).^n.*h0;
    f0 = h0;
    f1 = -h1;

    vk0 = conv(upsample(v00, L), f0); 
    vk1 = conv(upsample(v01, L), f1);
    
    v0 = sum([vk0; vk1]); 
end

function [y_n, h0, h1, f0, f1, H0, H1, F0, F1] = QMF(M, x_n)
    
    

    y1 = filter(f0, 1, upsample( dowsample( filter(h0, 1, x_n), M ), M) );
    y2 = filter(f1, 1, upsample( dowsample( filter(h1, 1, x_n), M ), M) );
        
    yk = [y1; y2];

    y_n = sum(yk);
end

function v = delay(x, n_D)

    v = [zeros(1, n_D), x];
end