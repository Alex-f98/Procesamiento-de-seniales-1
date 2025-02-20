% a y b)
clc; clear; close all;

N = 500;
sigmaW2 = 50;
wk = sqrt(sigmaW2)*randn(1,N);
k = 1:N;

a = 0.05;
b = 30 + 4*cos(2*pi*0.003*k);
fk = a*k + b;  %a y b constantes
yk = fk + wk;
%estas ecuaciones conforman la observacion yk = a.k + b + wk
%los estados salen de que a y b son constantes: X_k+1 = I.Xk
%x1 = a ,   x2 = b, SS?

F = eye(2); G = 0;
H = [1 1];%[k 1];
R = sigmaW2;
P0 = diag([20 500]);
X0 = zeros(2, 1);

%G solo por completitud
[x_est, errors] = Kalman( yk, P0, 0, X0, F, G , H, R);


figure()
    hold on
    plot(fk,'-o', 'LineWidth', 1, 'DisplayName', 'fk = ak + b')
    plot(x_est(1,:).*k + x_est(2,:) ,'-o', 'LineWidth', 1, 'DisplayName', 'fk_ = a_.k + b_')%'\hat{f}_k = \hat{a}.k + \hat{b}')
    
    plot(yk, 'DisplayName', 'señal observada: y_k')
     
    grid on
    title("Estimacion con Kalman")
    xlabel("Nro. iteraciones")
    ylabel("estimaciones")
    legend show
    
figure()
    r = xcorr(errors);
    [~, LL] = size(r);
    plot(r)
    title("Autocorrelacion de las innovaciones")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")

%% c)
%G solo por completitud
%eps_b = 0.5*randn(1,N)
G = eye(2);
Q = [0  0  ;
     0 0.5];
[x_est, errors] = Kalman( yk, P0, Q, X0, F, G , H, R);


figure()
    hold on
    plot(fk,'-o', 'LineWidth', 1, 'DisplayName', 'fk = ak + b')
    plot(x_est(1,:).*k + x_est(2,:) ,'-o', 'LineWidth', 1, 'DisplayName', 'fk_ = a_.k + b_')%'\hat{f}_k = \hat{a}.k + \hat{b}')
    
    plot(yk, 'DisplayName', 'señal observada: y_k')
     
    grid on
    title("Estimacion con Kalman")
    xlabel("Nro. iteraciones")
    ylabel("estimaciones")
    legend show
    
figure()  %esta maaaal!
    r = xcorr(errors);
    [~, LL] = size(r);
    plot( r)
    title("Autocorrelacion de las innovaciones")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")

    


    
    
    
    
%%
function [x_est, errors] = Kalman( Y, P0, Qk, X0, Fk, Gk ,Hk, Rk)

    N = length(Y);
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

    %aca poner asssert's para verificar que los tamaños sean validos.
    [~, n] = size(Fk);  [~, m] = size(Gk);   [l, ~] = size(Hk);
    


    Xk = X0;               %x_0/-1 = E[x0]
    Pk = P0;               %P_0/-1 = Cov[x - x0]

    errors = zeros(l, N);
    x_est = zeros(n, N);
    
    
    H = @(k)([k 1]);       %<--- parche para H variante en el tiempo.
    
    for k = 1:N

        %Actualizacion:
        K = Pk*H(k)'/(H(k)*Pk*H(k)' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/-1*H'_k + R_k)
        ek = Y(k) - H(k)*Xk;            %ek = yk - H.X_k/k-1)
        Xk = Xk + K*ek;                 %Tiene que ser un proceso Blanco.
        Pk = (eye(n,n) - K*H(k))*Pk;    %P_k/k = (I - K_k*H_k)*P_k/k-1 

        %Almaceno datos:
        errors(:, k) = ek;
        x_est(:, k) = Xk;


        %Prediccion:
        Xk = Fk*Xk;                 %X_k+1/k = F_k*X_k/k
        Pk = Fk*Pk*Fk'+ Gk*Qk*Gk'; %P_k+1/k = F_k*P_k/k*F'_k + G_k.Q_k.G'_k


    end
end

