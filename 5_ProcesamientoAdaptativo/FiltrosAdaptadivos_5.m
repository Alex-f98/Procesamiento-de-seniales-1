% RLS se quiere resolver el problema de identificasion del ej 3) con RLS
clc; clear; close all;
%Se quiere identificar un sistema de comunicaciones a través del cial se
%transmiten datos y(n) que atraviezan un canal con respuesta impulsiva
%h(n)={1 0.4 0.3 0.1 -0.2 0.05} mas ruido blanco gaussiano de media cero y
%varfianza sigmav2=0.02. la salida del sistema x= h*y|_N + v
h= [1 0.4 0.3 0.1 -0.2 0.05];
M = length(h);
sigmav2 = 0.02;

%) y(n) = 2*b(n) -1, b(n)~bernoulli(0.5), L = 2000. 
L = 2000;                       p = 0.5;
bn = binornd(1, p, L, 1);       vn = sqrt(sigmav2)*randn(L,1);

yn = 2*bn - 1;
% el filtro es: y(n)->|h(n)+v(n)|-->x(n)
y_h = filter(h, 1, yn);
xn = y_h + vn;

% a) Implemente el algoritmo RLS para la determinacion de los coeficientes
% Wn. Considere lambda = 0.99, delta = 0.001. Grafique los coeficientes del
% filtro estimado  Wn_ en funcion  de las iteraciones
delta = 0.001;
lambda = 0.99;

[w_, errors, Weigths] = RLS(yn, xn, M, delta, lambda);

figure()
stem(y_h, 'LineWidth', 2)
xlim([0,100])
title("x(n) vs n")

%Grafico los coeficientes
figure()
for k = 1:M
    hold on
    plot(Weigths(k,:), 'LineWidth', 2)
    plot([0, L], [h(k) h(k)], '-k')
end
grid on
ylabel("Estimacion de canal h(n):W(n)")
xlabel("Nro de iteracion")
title("Coeficientes h=w tal que N=M, entonces Jmin=sigmav2")

lgd = legend('$w_1$','','$w_2$','','$w_3$','','$w_4$','','$w_5$','','$w_6$','');
set(lgd,'Interpreter','latex', 'FontSize', 12);

%% b) repita la simulacion para al menos 500 realizaciones, calcule y grafique 
% la curva de aprendizaje Jn_ = 1/m sum|(eps_n)|^2

%Grafico la curva de aprendizaje
MM = 500;
JJ = zeros(MM, L-M+1);
for i = 1:MM
    
    bn = binornd(1, p, L, 1);
    yn = 2*bn - 1;

    vn = sqrt(sigmav2)*randn(L,1);

    % el filtroes: y(n)->|h(n)+v(n)|-->x(n)
    y_h = filter(h, 1, yn);
    xn = y_h + vn;
    
   [~, errors, Weigths] = RLS(yn, xn, M, delta, lambda);
   JJ(i,:) = abs(errors).^2 ;
end
%Curva de aprendizaje.
J = mean(JJ, 1);

figure()
semilogy(1:L-M+1, J)
title("Curva de aprendizaje con 500 iteraciones")
xlabel("m [ #iteraciones ]")
ylabel("$\hat{J}(n)$",'Interpreter','latex')
grid on

%% d) Para la misma señal x(n) del punto anterior, encuentre  los coeficientes 
% del filtro mediante LMS y RLS y compare las curvas de aprendizaje.
% considere lambda = 0.995 para RLS  y mu = 0.1 para LMS.
mu_LMS = 0.1;                   %
lambdaRLS = 0.995;              %factor de olvido
JJ_LMS = zeros(MM, L-5);
JJ_RLS = zeros(MM, L-5);
n = 0:L-1;

for i = 1:MM
    
    bn = binornd(1, p, L, 1);   vn = sqrt(sigmav2)*randn(L,1);
    yn = 2*bn - 1;
    % el filtroes: y(n)->|h(n)+v(n)|-->x(n)
    xn = filter(h, 1, yn) + vn;
    
   [WeigthsLMS, errorsLMS, ~] = LMS(yn, xn, M, mu_LMS);
   [~, errorsRLS, WeigthsRLS] = RLS(yn, xn, M, delta, lambdaRLS);
   
   JJ_LMS(i,:)  = abs(errorsLMS).^2 ;
   JJ_RLS(i,:) = abs(errorsRLS).^2 ;
end

%Curva de aprendizaje.
J_LMS = mean(JJ_LMS, 1);
J_RLS = mean(JJ_RLS, 1);

figure()
    hold on
    plot(1:L-5, log(J_LMS), 'b')
    plot(1:L-5, log(J_RLS), 'r')
    title("Curva de aprendizaje con 500 iteraciones $\mu = 0.1$, $\lambda = 0.995$", 'Interpreter', 'latex')
    xlabel("m [ #iteraciones ]")
    ylabel("$\hat{J}(n)$",'Interpreter','latex')
    legend("LMS", "RLS")
    grid on
    
 %RLS parece tener una curva de aprendizaje mas rapida y con menor error,
 %debido a que no es estocastico, sino deterministico.
%% e) repita el punto anterior con lambda = 0.95. Analice diferencias.
mu_LMS = 0.1;                   %
lambdaRLS = 0.95;              %factor de olvido
JJ_LMS = zeros(MM, L-5);
JJ_RLS = zeros(MM, L-5);

for i = 1:MM
    
    bn = binornd(1, p, L, 1);   vn = sqrt(sigmav2)*randn(L,1);
    yn = 2*bn - 1;
    % el filtroes: y(n)->|h(n)+v(n)|-->x(n)
    xn = filter(h, 1, yn) + vn;
    
   [WeigthsLMS, errorsLMS, ~] = LMS(yn, xn, M, mu_LMS);
   [~, errorsRLS, WeigthsRLS] = RLS(yn, xn, M, delta, lambdaRLS);
   
   JJ_LMS(i,:)  = abs(errorsLMS).^2 ;
   JJ_RLS(i,:) = abs(errorsRLS).^2 ;
end

%Curva de aprendizaje.
J_LMS = mean(JJ_LMS, 1);
J_RLS = mean(JJ_RLS, 1);

figure()
    hold on
    plot(1:L-5, log(J_LMS), 'b')
    plot(1:L-5, log(J_RLS), 'r')
    title("Curva de aprendizaje con 500 iteraciones $\mu = 0.1$, $\lambda = 0.95$", 'Interpreter', 'latex')
    xlabel("m [ #iteraciones ]")
    ylabel("$\hat{J}(n)$",'Interpreter','latex')
    legend("LMS", "RLS")
    grid on
%En este caso solo cambia lambda, es menos al anterior.
% esto quiere decir que a menor factor de olvido mayor sera lo que recuerde
% P de los pesos anteriores, haciendo que los primeros pesos influyan en el
% resulta final, aumentando en nivel de error que quedo, sabemos que el
% menor error va a ser el de la señal de ruido sigmav2, ~3.9db y a eso se
% le tienen que sumar los ruidos estocasticos del LMS que son mayores o los
% del algortimo RLS que son menores, pero que aumentan a medida que lambda
% se acerca mas a 1, es decir memoria infinica, lo importante es darle
% prioridad a las nuevas muestras.

%% g) perturbacion de h(n) el filtro en la mitad, graficar con distintos lambda
h1= [1 0.4 0.3 0.1 -0.2 0.05]; h2= [1 0.2 -0.5 0.3 0.02 0.1]; 
M = length(h);
% y1 = filter(h1, 1, yn); y2 = filter(h2, 1, yn);
% yh = [y1(1:L/2); y2(L/2+1:end)];
% xn = yh + vn;

lambdas = [0.5 0.9 0.995 1];
MM = 500;
JJ = zeros(MM, L-M+1);
    
for lambda = lambdas
   %Grafico la curva de aprendizaje
    
    for i = 1:MM

        bn = binornd(1, p, L, 1);       vn = sqrt(sigmav2)*randn(L,1);
        yn = 2*bn - 1;
        % el filtroes: y(n)->|h(n)+v(n)|-->x(n)
        y1 = filter(h1, 1, yn); y2 = filter(h2, 1, yn);
        y_h = [y1(1:L/2); y2(L/2+1:end)];
        xn = y_h + vn;

       [~, errors,~] = RLS(yn, xn, M, delta, lambda);
       JJ(i,:) = abs(errors).^2 ;
    end
    %Curva de aprendizaje.
    J = mean(JJ, 1);

    figure(10)
    hold all
    plot(1:L-M+1, log(J), 'DisplayName', "\lambda = " + lambda)
    title("Curva de aprendizaje con RLS 500 iteraciones")
    xlabel("m [ #iteraciones ]")
    ylabel("$\hat{J}(n)$",'Interpreter','latex')
    grid on 

end
    legend show
    
%Notar que a mayor lambda menor es el error, pero cuando hay una
%perturbacion el de mayor lambda tarda mas en volver al minimo error, puesto que tiene mayor memoria.

%% f) reinicie la matriz Pn en el instante justo en que se genera la perturbacion 
% del sistema. Observe las respuestas y saque conclusiones.
    

h1= [1 0.4 0.3 0.1 -0.2 0.05]; h2= [1 0.2 -0.5 0.3 0.02 0.1]; 
M = length(h);
% y1 = filter(h1, 1, yn); y2 = filter(h2, 1, yn);
% yh = [y1(1:L/2); y2(L/2+1:end)];
% xn = yh + vn;

lambdas = [0.5 0.9 0.995 1];
MM = 500;
JJ = zeros(MM, L-M+1);
    
for lambda = lambdas
   %Grafico la curva de aprendizaje
    
    for i = 1:MM

        bn = binornd(1, p, L, 1);       vn = sqrt(sigmav2)*randn(L,1);
        yn = 2*bn - 1;
        % el filtroes: y(n)->|h(n)+v(n)|-->x(n)
        y1 = filter(h1, 1, yn); y2 = filter(h2, 1, yn);
        y_h = [y1(1:L/2); y2(L/2+1:end)];
        xn = y_h + vn;

       [~, errors,~] = RLS2(yn, xn, M, delta, lambda, L/2);
       JJ(i,:) = abs(errors).^2 ;
    end
    %Curva de aprendizaje.
    J = mean(JJ, 1);

    figure(10)
    hold all
    plot(1:L-M+1, log(J), 'DisplayName', "\lambda = " + lambda)
    title("Curva de aprendizaje con RLS 500 iteraciones")
    xlabel("m [ #iteraciones ]")
    ylabel("$\hat{J}(n)$",'Interpreter','latex')
    grid on 

end
    legend show
    
    
    
    
    
    
    
    
    
    
    
    
    
%% funcion RLS

function [w_, errors, Weigths] = RLS(y, xn, M, delta, lambda)
L = length(y);
errors = zeros(1,L-M+1);
Weigths = zeros(M, L-M+1);
wn_ = zeros(M,1);                                        %w0_  0
Pn = eye(M)/delta;                                       %P0 = delta^-1.I_MxM
%Resto de las iteraciones: n: 1, 2,..., L-1.
for n = 1:L-M
   yn = flip(y(n:n+M-1));
   
   Kn = lambda\Pn*yn/(1 + lambda\yn'*Pn*yn);             %Ganancia(Mx1)
   %Kn = (Pn*yn/lambda)/(1 + yn'*Pn*yn/lambda);             %Ganancia(Mx1)
   eps = xn(n+M-1) - wn_'*yn;                            %Error a priori
   wn_ = wn_ + Kn * conj(eps);                           %Actualizacion
   Pn = lambda\Pn - lambda\Kn*yn'*Pn;                    %Mat inv de correlacion
   %Pn = Pn/lambda - Kn*yn'*Pn/lambda;                    %Mat inv de correlacion
   
   errors(n) = eps;
   Weigths(:, n+1) = wn_; %el primer peso estará asignado.
   
end
    w_ = wn_; 
end




%Funcion LMS, no es NLMS.
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









function [w_, errors, Weigths] = RLS2(y, xn, M, delta, lambda, LIMIT)
L = length(y);
errors = zeros(1,L-M+1);
Weigths = zeros(M, L-M+1);
wn_ = zeros(M,1);                                        %w0_  0
Pn = eye(M)/delta;                                       %P0 = delta^-1.I_MxM
%Resto de las iteraciones: n: 1, 2,..., L-1.
for n = 1:L-M
   yn = flip(y(n:n+M-1));
   
   Kn = lambda\Pn*yn/(1 + lambda\yn'*Pn*yn);             %Ganancia(Mx1)
   %Kn = (Pn*yn/lambda)/(1 + yn'*Pn*yn/lambda);             %Ganancia(Mx1)
   eps = xn(n+M-1) - wn_'*yn;                            %Error a priori
   wn_ = wn_ + Kn * conj(eps);                           %Actualizacion
   Pn = lambda\Pn - lambda\Kn*yn'*Pn;                    %Mat inv de correlacion
   %Pn = Pn/lambda - Kn*yn'*Pn/lambda;                    %Mat inv de correlacion
   if floor(LIMIT) == n
      Pn = eye(M)/delta;  
   end
   errors(n) = eps;
   Weigths(:, n+1) = wn_; %el primer peso estará asignado.
   
end
    w_ = wn_; 
end



