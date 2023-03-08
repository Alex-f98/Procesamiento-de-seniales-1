%Se quiere identificar un sistema de comunicaciones a través del cial se
%transmiten datos y(n) que atraviezan un canal con respuesta impulsiva
%h(n)={1 0.4 0.3 0.1 -0.2 0.05} mas ruido blanco gaussiano de media cero y
%varfianza sigmav2=0.02. la salida del sistema x= h*y|_N + v
clc; clear; close all;
h= [1 0.4 0.3 0.1 -0.2 0.05];
sigmav2 = 0.02;

%a) y(n) = 2*b(n) -1, b(n)~bernoulli(0.5), L = 2000. 
L = 2000;
p = 0.5;

bn = binornd(1, p, L, 1);
yn = 2*bn - 1;

vn = sqrt(sigmav2)*randn(L,1);

% el filtroes: y(n)->|h(n)+v(n)|-->x(n)
y_h = filter(h, 1, yn);
xn = y_h + vn;

figure()
stem(y_h, 'LineWidth', 2)
xlim([0,100])
title("x(n) vs n")
%% b) LMS con mu = 0.1, estima W
mu = 0.1;
N = length(h); %largo del filtro.
M = N;
[Weigths, errors, x_est] = LMS(yn, xn, M, mu);

%Grafico los coeficientes
figure()
for k=1:M
    hold on
    plot(Weigths(k,:), 'LineWidth', 2)
    plot([0, L], [h(k) h(k)], '-k')
end
grid on
title("Coeficientes h=w tal que N=M, entonces Jmin=sigmav2")

lgd = legend('$w_1$','','$w_2$','','$w_3$','','$w_4$','','$w_5$','','$w_6$','');
set(lgd,'Interpreter','latex', 'FontSize', 12);
%% c) Curva de aprendizaje 100 realizaciones
%Grafico la curva de aprendizaje
MM = 100;
JJ = zeros(MM, L-5);
for i = 1:MM
    
    bn = binornd(1, p, L, 1);
    yn = 2*bn - 1;

    vn = sqrt(sigmav2)*randn(L,1);

    % el filtroes: y(n)->|h(n)+v(n)|-->x(n)
    y_h = filter(h, 1, yn);
    xn = y_h + vn;
    
   [Weigths, errors, x_est] = LMS(yn, xn, M, mu);
   JJ(i,:) = abs(errors).^2 ;
end
%Curva de aprendizaje.
J = mean(JJ, 1);

figure()
semilogy(1:L-5, J)
title("Curva de aprendizaje con 100 iteraciones")
xlabel("m [ #iteraciones ]")
ylabel("$\hat{J}(n)$",'Interpreter','latex')
grid on

%% d) Curva de aprendizaje 100 realizaciones u = 0.001
%Grafico la curva de aprendizaje
mu = 0.005;
MM = 100;
JJ = zeros(MM, L-5);
for i = 1:MM
    
    bn = binornd(1, p, L, 1);
    yn = 2*bn - 1;

    vn = sqrt(sigmav2)*randn(L,1);

    % el filtroes: y(n)->|h(n)+v(n)|-->x(n)
    y_h = filter(h, 1, yn);
    xn = y_h + vn;
    
   [Weigths, errors2, x_est] = LMS(yn, xn, M, mu);
   JJ(i,:) = abs(errors2).^2 ;
end
%Curva de aprendizaje.
J = mean(JJ, 1);

figure()
semilogy(1:L-5, J)
title("Curva de aprendizaje con 100 iteraciones mu = 0.001")
xlabel("m [ #iteraciones ]")
ylabel("$\hat{J}(n)$",'Interpreter','latex')
grid on


%%
function [Weigths, errors, x_est] = LMS(yn, x, N, u, w_inicial) 
%yn(vector columna)
%w_inicial (vector columna)
assert( (size(yn,2)==1 || size(yn,2)==1),"x(n) e y(n) tienen que ser vectores columna");

if ~exist('w_inicial')
        w_inicial = zeros(N,1); 
end
if size(w_inicial,2)~=1; w_inicial = w_inicial'; end
Weigths =[];  
errors = [];
x_est =[];

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
   for i = 1:L-N+1
        Y_    = flip(yn(i:i+N-1));
        xn_   = w' * Y_;
        error = x(i+N-1) - xn_;
        w     = w + u * Y_ * error; %conj(error);
        
        Weigths(:, end+1) = w;
        errors(end+1) = error;
        x_est(end+1) = xn_; 
   end
   
end


