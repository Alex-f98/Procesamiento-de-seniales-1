%1) estimacion de una constante
clc; clear; close all;
% b se mantiene constante en el tiempo. b_k+1 = b_k, se quiere estimar el
% valor a partir de medicions ruidosas yk = bk + wk / wk ~ N(0; 0.25).

%a) encontrar Fd y H.

%ya está discretizado.
%X' = F X + G u         y = H x
Fd = 1;  G = 0;
H = 1;

%b) Genere un proceso de mediciones yk con b=2 N=500;
N = 500;
bk = 2;
sigma_W2 = 0.25;

wk = sqrt(sigma_W2)*randn(N, 1);
yk = H*bk + wk;

%% c) implementar Kalman (guardar las innovaciones ek = yk - H.X_k/k-1)
% x_0/-1 = 0; P_0/-1 = 10; 

% n: cantidad de estados xk
% m: dimensión del ruido del proceso vk
% l: cantidad de ecuaciones de medición yk
% x: Estados (nx1)
% y: Ecuaciones de medición (lx1)
% v: Ruido del proceso (mx1)
% w: Ruido de medición (lx1)

% F (nxn)             G (nxm)
% H (lxn)             R (lxl)
% Q (mxm)     P (nxn)         K (nxl)
Y = yk;
po = 10;
Fk = Fd;
Gk = G;
Hk = H;
Rk = sigma_W2;

[n, ~] = size(Hk);
[~, m] = size(G);
[l, ~] = size(Hk);

Xk = zeros(n, 1);      %x_0/-1 = 0
Pk = po*ones(n, n);    %P_0/-1 = 10;

errors = zeros(l, N);
x_est = zeros(n, 1);
for k = 1:N
    
    %Actualizacion:
    K = Pk*Hk'/(Hk*Pk*Hk' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/-1*H'_k + R_k)
    ek = Y(k) - Hk*Xk;            %ek = yk - H.X_k/k-1)
    Xk = Xk + K*ek;               %Tiene que ser un proceso Blanco.
    Pk = (eye(n,n) - K*Hk)*Pk;    %P_k/k = (I - K_k*H_k)*P_k/k-1 
    
    %Almaceno datos:
    errors(:, k) = ek;
    x_est(:, k) = Xk;
    
    
    %Prediccion:
    Xk = Fk*Xk;                 %X_k+1/k = F_k*X_k/k
    Pk = Fk*Pk*Fk'; %+ Gk*Qk*Gk'; %P_k+1/k = F_k*P_k/k*F'_k + G_k.Q_k.G'_k
    
    
end

figure()
    hold on
    plot(x_est, '-k', 'LineWidth', 2, 'DisplayName', 'Señal estimada: X_{est}')
    plot(yk, '-b', 'DisplayName', 'señal observada: y_k')
    plot([1 N], [bk bk], '-r', 'LineWidth', 1, 'DisplayName', 'Valor real: b')
    grid on
    title("Estimacion de b = 2, con Kalman")
    xlabel("Nro. iteraciones")
    ylabel("Constante b")
    legend show
    
figure()
    r = xcorr(errors);
    plot(-floor(length(r)/2):floor(length(r)/2) , r)
    title("Autocorrelacion de las innovaciones, PROCESO BLANCO")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")

