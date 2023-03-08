%a)ALGORITMO LMS.
clc; clear;close all;
%b)
L = 2000;
Sigmav2 = 0.1;
M = 2;
u = 0.01;

h_y = [1 0.5 0.25];
h_s = [1 1.2 0.60 0.3];  %<--- habia un 0 en vez de 1.

fn = randn(1, L);       %<---rand() en ves de randn()
vn = sqrt(Sigmav2)*randn(1, L);
sn = filter(h_s, 1, fn);

yn = filter(h_y, 1, fn);
xn = sn + vn;

[Weigth, errors] = LMS(yn.', xn.', M, u);

%Grafico los coeficientes
figure()
hold on
plot(Weigth(1,:), 'r', 'LineWidth', 2)
plot(Weigth(2,:), 'b', 'LineWidth', 2)
grid on

%Grafico la curva de aprendizaje
MM = 500;
JJ = zeros(MM, L-1);
for i = 1:MM
    
    fn = randn(1, L);       %<---rand() en ves de randn()
    vn = sqrt(Sigmav2)*randn(1, L);
    sn = filter(h_s, 1, fn);

    yn = filter(h_y, 1, fn);
    xn = sn + vn;
    
    [~, errors] = LMS(yn.', xn.', M, u) ;
   JJ(i,:) = abs(errors).^2 ;
end
%Curva de aprendizaje.
J = mean(JJ, 1);

figure()
semilogy(1:L-1, J)
title("Curva de aprendizaje con 200 iteraciones")
xlabel("m [ #iteraciones ]")
ylabel("$\hat{J}(n)$",'Interpreter','latex')
grid on

%% c)
% % Curvas de nivel de J(w0, w1)
w0    = [-3.19; 4.47];
w1    = [ 4.84; 4.52];
w2    = [ 1.65; 8.99];
wcond = [w0, w1, w2];
sigma2_x = var(xn);
nW = 9;            %cantidad de puntos de la grilla.

yy = zeros(M,M);
yx = zeros(M,1);
% 1. Estimar R utilizando ventanas de largo M > K
for j = 1: L-M+1
   yM = transpose(flip(yn(j:j+M-1)));

   yy = yy + yM*yM';
   yx = yx + yM*xn(j+M-1);
end
%estimacion de Ry
Ry_  = yy/M;                      %matriz de correlacion MxM
%Estimacion de Ryx
Ryx_ = yx/M;

[W1, W2] = meshgrid(linspace(-nW, nW), linspace(-nW, nW));

J = zeros(size(W1));
%Recordar: J(W) = sigma2x + W'Ryx - Ryx'W + W'RyW
for idx = 1:length(W1)^2
    W = [W1(idx); W2(idx)];
    J(idx) = sigma2_x + W'*Ryx_ - Ryx_'*W + W'*Ry_*W;
end


MM = 500;
JJ = zeros(MM, L-1);

for cond = 1:3
    for i = 1:MM

        fn = randn(1, L);       %<---rand() en ves de randn()
        vn = sqrt(Sigmav2)*randn(1, L);
        sn = filter(h_s, 1, fn);

        yn = filter(h_y, 1, fn);
        xn = sn + vn;

        [Weigths, errors] = LMS(yn.', xn.', M, u, wcond(:, cond)) ;
    end
    
    figure(7) 
    contour(W1, W2, J, 50)
    hold on
    plot(Weigths(1,:), Weigths(2,:), '-*');
    
end
    legend("", "$w_0$","", "$w_1$","", "$w_2$",'Interpreter','latex')
    xlabel('$w_1$','interpreter','latex','FontSize',12);
    ylabel('$w_2$','interpreter','latex','FontSize',12);
    title("Trayectoria de los pesos para distintos pesos inciales")
%% d)

hold off

mu = 0.2;

% % Curvas de nivel de J(w0, w1)
w0    = [-3.19; 4.47];
w1    = [ 4.84; 4.52];
w2    = [ 1.65; 8.99];
wcond = [w0, w1, w2];
sigma2_x = var(xn);
nW = 9;            %cantidad de puntos de la grilla.

yy = zeros(M,M);
yx = zeros(M,1);
% 1. Estimar R utilizando ventanas de largo M > K
for j = 1: L-M+1
   yM = transpose(flip(yn(j:j+M-1)));

   yy = yy + yM*yM';
   yx = yx + yM*xn(j+M-1);
end
%estimacion de Ry
Ry_  = yy/M;                      %matriz de correlacion MxM
%Estimacion de Ryx
Ryx_ = yx/M;

[W1, W2] = meshgrid(linspace(-nW, nW), linspace(-nW, nW));

J = zeros(size(W1));
%Recordar: J(W) = sigma2x + W'Ryx - Ryx'W + W'RyW
for idx = 1:length(W1)^2
    W = [W1(idx); W2(idx)];
    J(idx) = sigma2_x + W'*Ryx_ - Ryx_'*W + W'*Ry_*W;
end


MM = 500;
JJ = zeros(MM, L-1);

for cond = 1:3
    for i = 1:MM

        fn = randn(1, L);       %<---rand() en ves de randn()
        vn = sqrt(Sigmav2)*randn(1, L);
        sn = filter(h_s, 1, fn);

        yn = filter(h_y, 1, fn);
        xn = sn + vn;

        [Weigths, errors] = LMS(yn.', xn.', M, mu, wcond(:, cond)) ;
    end
    
    figure(8) 
    contour(W1, W2, J, 50)
    hold on
    plot(Weigths(1,:), Weigths(2,:), '-*');
    
end
    legend("", "$w_0$","", "$w_1$","", "$w_2$",'Interpreter','latex')
    xlabel('$w_1$','interpreter','latex','FontSize',12);
    ylabel('$w_2$','interpreter','latex','FontSize',12);
    title("Trayectoria de los pesos para distintos pesos inciales, $\mu = 0.2$",'interpreter','latex')

%%
function [Weigths, errors] = LMS(yn, x, N, u, w_inicial) 
%yn(vector columna)
%w_inicial (vector columna)
assert( (size(yn,2)==1 || size(yn,2)==1),"x(n) e y(n) tienen que ser vectores columna");

if ~exist('w_inicial')
        w_inicial = zeros(N,1); 
end
if size(w_inicial,2)~=1; w_inicial = w_inicial'; end
Weigths=[];  
errors = [];
%XPICO=[];

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
        %XPICO(end+1) = xn_; 
   end
   
end
