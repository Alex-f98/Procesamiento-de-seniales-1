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
sigma_b2 = 0.2;                 %A MAYOR Q, CONFIO MAS EN LAS MEDICIONES QUE EN EL MODELO, POR ESO ES MAS RUIDOSO.
%Fd = e^F.T;
Fd = eye(2) + F*T + (F*T)^2/2;
Gd = eye(2,2);                  %G (nxm)

%Qd = E[Vk.Vk^H]:
Q = [sigma_b2 0; 0 0];          %Cov(eps_v * eps_v')
Qd = Q*T;                       % Q (mxm)

% yk = x2 + wk
Hd = [0 1];
%Rd = sigmaW2;
wk = sqrt(sigmaW2)*randn(1, length(v));

%a)
yk = Hd*[v, d]' + wk;

P0 = diag([20 1]);
X0 = zeros(2, 1);
%b)Rd= 1,2;

Rd = 1.2;
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
    
%% c)
Rd = 0.00025;
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
    
%% d)

for sigmaW2 = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5]
    Rd = sigmaW2;
    
    [~, ~, Ks] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);

 figure(1)
    subplot(1, 2, 1)
    hold on
    plot(Ks(1, :), 'LineWidth', 2, 'DisplayName', "\sigma_W^2 = "+sigmaW2);
    hold off
    subplot(1, 2, 2)
    hold on
    plot(Ks(2, :), 'LineWidth', 2, 'DisplayName', "\sigma_W^2 = "+sigmaW2);
    hold off
end
    
    subplot(1,2,1)
    xlabel("Nro iteraciones")
    ylabel("K1")
    grid minor
    legend show
    ylim([-0.01, 0.7])
    subplot(1,2,2)
    xlabel("Nro iteraciones")
    ylabel("K2")
    grid minor
    legend show
    ylim([-0.01, 0.17])
    title("Ganancias del FK variando Q, mayor Q mayor K")
    
