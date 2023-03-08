clc; clear; close all;

L = 2000;
sigmav2 = 0.02;
h= [1 0.4 0.3 0.1 -0.2 0.05];

mu_LMS = 0.005;
mu_NLMS = 0.4;

MM = 500;
%%
N = length(h); %largo del filtro.
M = N;
JJ_LMS = zeros(MM, L-5);
JJ_NLMS = zeros(MM, L-5);

for i = 1:MM
    
    yn = sqrt(36)*randn(L,1);
    vn = sqrt(sigmav2)*randn(L,1);

    % el filtro es: y(n)->|h(n)+v(n)|-->x(n)
    y_h = filter(h, 1, yn);
    xn = y_h + vn;
    
   [WeigthsLMS, errorsLMS, ~] = LMS(yn, xn, M, mu_LMS);
   [WeigthsNLMS, errorsNLMS, ~] = NLMS(yn, xn, M, mu_NLMS);
   
   JJ_LMS(i,:)  = abs(errorsLMS).^2 ;
   JJ_NLMS(i,:) = abs(errorsNLMS).^2 ;
end
%Curva de aprendizaje.
J_LMS = mean(JJ_LMS, 1);
J_NLMS = mean(JJ_NLMS, 1);

figure(1)
    hold on
    plot(1:L-5, log(J_LMS))
    plot(1:L-5, log(J_NLMS))
    title("Curva de aprendizaje con 100 iteraciones mu = 0.001")
    xlabel("m [ #iteraciones ]")
    ylabel("$\hat{J}(n)$",'Interpreter','latex')
    legend("LMS", "NLMS")
    grid on
    
%% c) 1/|y(n)|^2 se peude ir a infinito
clc; clear; close all;

L = 2000;
sigmav2 = 0.02;
h= [1 0.4 0.3 0.1 -0.2 0.05];
delta = 0; %0.001;

mu_LMS = 0.0025;%0.005;
mu_NLMS = 0.05;

MM = 500;

N = length(h); %largo del filtro.
M = N;
JJ_LMS = zeros(MM, L-5);
JJ_NLMS = zeros(MM, L-5);
n = 0:L-1;
yn = 1 + cos(0.004*pi*n)';

for i = 1:MM
    
    %sqrt(36)*randn(L,1);
    vn = sqrt(sigmav2)*randn(L,1);

    % el filtro es: y(n)->|h(n)+v(n)|-->x(n)
    y_h = filter(h, 1, yn);
    xn = y_h + vn;
    
   [WeigthsLMS, errorsLMS, ~] = LMS(yn, xn, M, mu_LMS);
   [WeigthsNLMS, errorsNLMS, ~] = NLMS(yn, xn, M, mu_NLMS);
   
   JJ_LMS(i,:)  = abs(errorsLMS).^2 ;
   JJ_NLMS(i,:) = abs(errorsNLMS).^2 ;
end
%Curva de aprendizaje.
J_LMS = mean(JJ_LMS, 1);
J_NLMS = mean(JJ_NLMS, 1);

figure(1)
    hold on
    plot(1:L-5, log(J_LMS))
    plot(1:L-5, log(J_NLMS))
    title("Curva de aprendizaje con 100 iteraciones mu = 0.001")
    xlabel("m [ #iteraciones ]")
    ylabel("$\hat{J}(n)$",'Interpreter','latex')
    legend("LMS", "NLMS")
    grid on
    

    



%% funciones
function [Weigths, errors, x_est] = NLMS(yn, x, N, u, w_inicial, delta) 
%yn(vector columna)
%w_inicial (vector columna)
assert( (size(yn,2)==1 || size(yn,2)==1),"x(n) e y(n) tienen que ser vectores columna");

if ~exist('w_inicial')
    w_inicial = zeros(N,1); 
end

% if ~exist('delta')
%    delta = 0;
% end
delta = 0.001;   %<----harcodeado!

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
        k = u/(delta + (Y_'*Y_));
        w     = w + k * Y_ * error; %conj(error);
        
        Weigths(:, end+1) = w;
        errors(end+1) = error;
        x_est(end+1) = xn_; 
   end
   
end




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

