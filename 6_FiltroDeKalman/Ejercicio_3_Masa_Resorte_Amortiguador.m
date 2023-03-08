%Ejercicio_3_Masa_resorte_amortiguador
clc; clear; close all;
load('MRA.mat')

%k_e = 100;  %N/m
%M = 10;     %Kg
%b = 10;     %Kg/s
%Aca es lo mismo que el 2) pero con b_=18 en ves de 10, hay ruido de
%proceso pero supuestamente no lo sé, error de modelado, Q = 0, no lo tengo
%en cuenta.
b_ = 18;     %Kg/s

%% a) asumiento que no hay fuerzas externas(u(t)=0): m.a(t) = -k_e.d(t) -b.v(t), hay cond. inicial
%x1 = v    x2=d (velocidad y posicion)

%[x1'  x2'] = [-b -k; 1 0][x1 x2]^T
%Esto es en tiempo continuo:
F = [-b_ /M -Ke/M; 1 0];
G = 0;

%% b) yk = dk + wk,  T = 10ms, wk~N(0; 0.05) es el ruido de midicion
sigmaW2 = 0.025;
%Fd = e^F.T;
Fd = eye(2) + F*T + (F*T)^2/2;      Gd = 0;
%Qd = E[Vk.Vk^H]:
Qd = 0;                                     %No asumimos ruido de proceso.

% yk = x2 + wk
Hd = [0 1];
Rd = sigmaW2;
wk = sigmaW2*randn(1, length(v));

yk = Hd*[v, d]' + wk;
P0 = diag([20 1]);
X0 = zeros(2, 1);


[x_est, errors] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);


figure()
    hold on
    plot(x_est(1,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(1)')
    plot(x_est(2,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(2)')
    
    plot(yk, 'DisplayName', 'señal observada: y_k')
    
    plot(d, 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = d(n)')
    plot(v, 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = v(n)')
    
    grid minor
    title("Estimacion con Kalman, b = 18 \neq 10, Q = 0")
    xlabel("Nro. iteraciones")
    ylabel("estimaciones")
    legend show


    
    
figure()
    r = xcorr(errors);
    [~, LL] = size(r);
    %plot(-floor(LL/2):floor(LL/2)-1 , r)
    plot(r)
    title("Autocorrelacion de las innovaciones")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")

%% b) ahora admitimos que hay un error en el modelo Qd = Q.T, Q es la covarianza de proceso continuo
% ajustar sigma_b2 hasta que hasta ajustar las variables de estado.

%ma = Fe + Fr, Fe = -ke.d,  Fr = -b.d'
%m.d(t)'' = -k_e.d(t) -b.v(t) + eps_v/ X' = F.X + G.U;  Y= H.X;

sigma_b2 = 1.5;                 %A MAYOR Q, CONFIO MAS EN LAS MEDICIONES QUE EN EL MODELO, POR ESO ES MAS RUIDOSO.
Gd = eye(2,2);                  %G (nxm)
Q = [sigma_b2 0; 0 0];          %Cov(eps_v * eps_v')
Qd = Q*T;                       % Q (mxm)
[x_est, errors] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);


figure()
    hold on
    plot(x_est(1,:),'-', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(1)')
    plot(x_est(2,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(2)')
    
    plot(yk, 'DisplayName', 'señal observada: y_k')
    
    plot(d, 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = d(n)')
    plot(v, 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = v(n)')
    
    grid minor
    title("Estimacion con Kalman, b = 18 \neq 10, \sigma_b ^2 = "+sigma_b2)
    xlabel("Nro. iteraciones")
    ylabel("estimaciones")
    legend show


    
    
figure()
    r = xcorr(errors);
    [~, LL] = size(r);
    %plot(-floor(LL/2):floor(LL/2)-1 , r)
    plot(r)
    title("Autocorrelacion de las innovaciones")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")

%% c)


%A MAYOR Q, CONFIO MAS EN LAS MEDICIONES QUE EN EL MODELO, POR ESO ES MAS RUIDOSO.

for sigma_b2 = [0.1 0.2 0.3 0.4 0.5]
    
    Q = [sigma_b2 0; 0 0];          %Cov(eps_v * eps_v')
    Qd = Q*T;                       % Q (mxm)
    
    [~, ~, Ks] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);
    
    figure(1)
    subplot(1, 2, 1)
    hold on
    plot(Ks(1, :), 'LineWidth', 2, 'DisplayName', "\sigma_b^2 = "+sigma_b2);
    hold off
    subplot(1, 2, 2)
    hold on
    plot(Ks(2, :), 'LineWidth', 2, 'DisplayName', "\sigma_b^2 = "+sigma_b2);
    hold off
    
end

    subplot(1,2,1)
    xlabel("Nro iteraciones")
    ylabel("K1")
    grid minor
    legend show
    subplot(1,2,2)
    xlabel("Nro iteraciones")
    ylabel("K2")
    grid minor
    legend show
    title("Ganancias del FK variando Q, mayor Q mayor K")
    











%%
function [x_est, innovaciones, Ks] = Kalman( Y, P0, Qk, X0, Fk, Gk ,Hk, Rk)

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

    [~, n] = size(Fk);  [~, m] = size(Gk);   [l, ~] = size(Hk);
    


    Xk = X0;               %x_0/-1 = E[x0]
    Pk = P0;               %P_0/-1 = Cov[x - x0]

    innovaciones = zeros(l, N);
    x_est = zeros(n, N);
    Ks = zeros(n, N);
   
    
    for k = 1:N

        %Actualizacion:
        K = Pk*Hk'/(Hk*Pk*Hk' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/k-1*H'_k + R_k)
        ek = Y(k) - Hk*Xk;            %ek = yk - H.X_k/k-1)
        Xk = Xk + K*ek;               %Tiene que ser un proceso Blanco.
        Pk = (eye(n,n) - K*Hk)*Pk;    %P_k/k = (I - K_k*H_k)*P_k/k-1 

        %Almaceno datos:
        innovaciones(:, k) = ek;
        x_est(:, k) = Xk;
        Ks(:, k) = K;


        %Prediccion:
        Xk = Fk*Xk;                 %X_k+1/k = F_k*X_k/k
        Pk = Fk*Pk*Fk' + Gk*Qk*Gk'; %P_k+1/k = F_k*P_k/k*F'_k + G_k.Q_k.G'_k


    end
end












