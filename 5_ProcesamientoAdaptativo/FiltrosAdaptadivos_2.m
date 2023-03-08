clc; clear; close all;
%a) el error optimo eo en este caso es la señal S(n) dado que con el filtro
%adaptativo estoy estimando x(n)_ ~ g(n)=Acos(wi.n + theta).

%b) g(n){ A~N(0,01; 0,001)     theta~U(-pi, pi)}(son VA pero toman un valor fijo en cada realizacion)
[sn, fs] = audioread('CANCIONES/Pista_05.wav');

%Generación de la interferencia
f0 = 1.5e3;   %[Hz] Frecuencia continua
fd = f0/fs; %Frecuencia discreta
wd = 2*pi*fd;
M = 2;
mu = 5e-5;

n = (1:length(sn))';
%Parametros estocasticos
A = normrnd(0.01,0.001);        % A ~ N(0.01;0.001)
theta = unifrnd(-pi,pi);        % theta ~ U(-pi,pi)
%Señal ruidosa
gn = A*cos(wd*n + theta);       % g(n)=A.cos(wi.n + theta)

xn = sn + gn;

yn = [cos(wd*n), sin(wd*n)];    %[ ||, || ]
%c) mu = 5e-5; graficar los M coeficientes del filtro en funcion de las iteraciones; las
%diferencias cuadraticas entre la señal de interferencia y salida del
%filtro LMS |g(n) - x(n)_|^2.

[w, Weigths, errors, X_est] = LMS2(yn.', xn, M, mu);

%Grafico los coeficientes
figure()
    hold on
    plot(Weigths(1,:), 'r', 'LineWidth', 2)
    plot(Weigths(2,:), 'b', 'LineWidth', 2)
    grid on

figure()
    plot((gn' - X_est).^2);
    grid on, grid minor;
    title('Diferencia cuadratica $|g(n)-\hat{x}(n)|$','Interpreter','latex');
    ylabel('$|g(n)-\hat{x}(n)|$','Interpreter','latex')
    xlabel('Nro de iteración','FontSize',12);

%% e) repetir c) con g2(n){g2 = g/ n<=L/2; g2=2*g/ n>L/2}
L = length(gn);
g2n = [gn(1:L/2); 2*gn(L/2+1:end)];   %piecewise()

x2n = sn + g2n;

[w, Weigths, errors, X_est] = LMS2(yn.', x2n, M, mu);

%Grafico los coeficientes
figure()
    hold on
    plot(Weigths(1,:), 'r', 'LineWidth', 2)
    plot(Weigths(2,:), 'b', 'LineWidth', 2)
    grid on

figure()
    plot((g2n' - X_est).^2);
    grid on, grid minor;
    title('Diferencia cuadratica $|g2(n)-\hat{x}(n)|$','Interpreter','latex');
    ylabel('$|g(n)-\hat{x}(n)|$','Interpreter','latex')
    xlabel('Nro de iteración','FontSize',12);

%% f) reproducir la señal contaminada y luego la filtrada.
sound(sn, fs)
%%
sound(x2n, fs)
%%
sound(errors, fs)
%se escucha como al inicio hay ruido, luego los pesos se van adaptado y se
%elimina el ruido, como este era el caso que cambiaba g2n en la mtad se
%puede escuchar como reaparece el ruido a la mitad y luego se filtra.
%% algortimo LMS para cancelador de ruido
% la diferencia es que mi señal de observada son senos y cosenos que
% modelan al ruido, asi mi filtro adaptativo acomoda sus pesos para que
% x_=g sea ruido luego mi señal de referencia sera s = s + g entonces el
% error sera la señal s, con ello se realimenta hasta que se cancele el
% ruido

function [w, Weigths, errors, X_est] = LMS2(yn, x, N, u, w_inicial)
%x: señal de referencia
%yn: señal obsservada
%N: largo de la ventana para estima Y_ con yn
%u: Constante mu del filtro.
%w_inicial: Inicializacion de los pesos del filtro(Opcional)

%yn(vector columna)
%w_inicial (vector columna)
%assert( (size(yn,2)==1 || size(yn,2)==1),"x(n) e y(n) tienen que ser vectores columna");

if ~exist('w_inicial')
        w_inicial = zeros(N,1); 
end
if size(w_inicial,2)~=1; w_inicial = w_inicial'; end
Weigths=[];  
errors = [];
X_est=[];

%u; %tiene un rango 0 < mu < 1/traza(A)
%a) ALGORITMO LMS
L = length(yn);
% ➢ definir las condiciones iniciales W(0)
    w = w_inicial;
% ➢ Paso 0: estimar W(1) a partir de W(0)
% ➢ Paso n: estimar W(n+1) a partir de W(n):
    % 1. Se calcula la salida del filtro,          x_(n)= W'_n.Y(n).
    % 2. Se calcula el error de estimación,        e(n) = x(n) - x_(n).
    % 3. Se adaptan los coeficientes del filtro,   W{n+1}_= Wn_ + mu.Y(n)e(n)*
   for i = 1:L%-N+1
        Y_    =  yn(:,i); %flip(yn(i:i+N-1));
        xn_   = w' * Y_;
        error = x(i) - xn_;
        w     = w + u * Y_ * error; %conj(error);
        
        Weigths(:, end+1) = w;
        errors(end+1) = error;
        X_est(end+1) = xn_; 
   end
   
end











