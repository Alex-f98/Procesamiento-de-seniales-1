
clc; clear; close all;
%Problema 1
load('tp5.mat'); 
[N, M] = size(p);
%     a: [351×3 double]
%     p: [351×3 double]
%     v: [351×3 double]

%p: posicion; v: velocidad;  a: aceleracion
T = 1;      %[s]

Fc = [0 0 0, 1 0 0, 0 0 0;
     0 0 0, 0 1 0, 0 0 0;
     0 0 0, 0 0 1, 0 0 0;%     
     0 0 0, 0 0 0, 1 0 0;
     0 0 0, 0 0 0, 0 1 0;
     0 0 0, 0 0 0, 0 0 1;%
     0 0 0, 0 0 0, 0 0 0;
     0 0 0, 0 0 0, 0 0 0;
     0 0 0, 0 0 0, 0 0 0 ];
 
q = diag([10^-2, 8e-3, 2e-5]);
      
Qc = zeros(size(Fc)); Qc(end-2:end, end-2:end) = q;

Fd = eye(size(Fc)) + Fc*T + (Fc*T)^2/2 ;      Gd = eye(size(Fc));
%Fd = eye(size(Fc)) + Fc*T + (Fc*T)^2/2 + (Fc*T)^3/3;      Gd = eye(size(Fc));

Qd = Qc*T;
%Qd = Qc*T + (Fc*Qc + 0.5*Qc*Fc')*T^2;

%Hd = eye(9,9); %[mxn]
Hd = zeros(3,9); Hd(:, 1:3) = eye(M);

Sigma_p2 = 2500;    %[m^2]
Sigma_v2 = 25;      %[m^2/s^2]
Sigma_a2 = 1;       %[m^2/s^4]

Rd = eye(M,M)*Sigma_p2;
%Rd = diag([ones(1,3)*Sigma_p2, ones(1,3)*Sigma_v2, ones(1,3)*Sigma_a2]); %[mxm]

nu_k = sqrt(Sigma_p2)*randn(M, N);

yk = Hd*[p, v, a]' + nu_k;

X0 = [10 -600 50 0 0 0 0 0 0]';
P0 = diag([10^6 10^6 10^6 10^4 10^4 10^4 5 5 5]);

%-%(a) Utilizando una matriz de covarianza de medici´on mayor que la real: R0 = 10^3.R.
Rd_a = Rd*10^3;
Rd_b = Rd*10^-3;
[x_est, innovaciones] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd_a);
[x_est_, innovaciones_] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd_b);

 figure(1)
subplot(2, 1, 1)
    hold on
    plot(x_est(1,:),'-o', 'LineWidth', 1, 'DisplayName', 'p_x: X_{est}(1)')
    plot(x_est(2,:),'-o', 'LineWidth', 1, 'DisplayName', 'p_y: X_{est}(2)')
    plot(x_est(3,:),'-o', 'LineWidth', 1, 'DisplayName', 'p_z: X_{est}(3)')
    
    plot(p(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = p_x')
    plot(p(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = p_y')
    plot(p(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_3 = p_z')
    
    plot(yk(1, :), 'LineWidth', 2, 'DisplayName', 'Valor medido: y_1')
    plot(yk(2, :), 'LineWidth', 2, 'DisplayName', 'Valor medido: y_2')
    plot(yk(3, :), 'LineWidth', 2, 'DisplayName', 'Valor medido: y_3')

    grid minor
    xlabel("Nro. iteraciones"); ylabel("estimaciones"); legend('show','FontSize', 13 )
    title("Estimacion con Kalman, R' = 10^3R")
    
    hold off
%(b) Utilizando una matriz de covarianza de medici´on menor que la real: R0 = 10−3R
subplot(2, 1, 2)
    hold on
    plot(x_est_(1,:),'-o', 'LineWidth', 1, 'DisplayName', 'p_x: X_{est}(1)')
    plot(x_est_(2,:),'-o', 'LineWidth', 1, 'DisplayName', 'p_y: X_{est}(2)')
    plot(x_est_(3,:),'-o', 'LineWidth', 1, 'DisplayName', 'p_z: X_{est}(3)')
    
    plot(p(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = p_x')
    plot(p(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = p_y')
    plot(p(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_3 = p_z')
    
    plot(yk(1, :), 'LineWidth', 2, 'DisplayName', 'Valor medido: y_1')
    plot(yk(2, :), 'LineWidth', 2, 'DisplayName', 'Valor medido: y_2')
    plot(yk(3, :), 'LineWidth', 2, 'DisplayName', 'Valor medido: y_3')

    grid minor
    xlabel("Nro. iteraciones"); ylabel("estimaciones"); legend('show','FontSize', 13 )
    title("Estimacion con Kalman, R' = 10^-3 R")
    hold off
    
%% innovaciones

figure()
    hold on;
    rx = xcorr(innovaciones(1,:));
    ry = xcorr(innovaciones(2,:));
    rz = xcorr(innovaciones(3,:));
    plot(rx, 'DisplayName', 'p_x'); plot(ry, 'DisplayName', 'p_y'); plot(rz, 'DisplayName', 'p_z')
    title("Autocorrelacion de las innovaciones , p R' = 10^3R")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")
    grid minor
    legend('show','FontSize', 15 )
    
figure()
    hold on;
    rx = xcorr(innovaciones_(1,:));
    ry = xcorr(innovaciones_(2,:));
    rz = xcorr(innovaciones_(3,:));
    plot(rx, 'DisplayName', 'p_x'); plot(ry, 'DisplayName', 'p_y'); plot(rz, 'DisplayName', 'p_z')
    title("Autocorrelacion de las innovaciones, p, R' = 10^-3 R ")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")
    grid minor
    legend('show','FontSize', 15 )
    

    
    
