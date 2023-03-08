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

%% (a) Para las posiciones, estados iniciales cercanos al real con una varianza inicial pequeña.
%x0/−1 = [−15 10 960 0 0 0 0 0 0]T
%P0/−1 = diag([10 10 10 104 104 104 5 5 5])

X0 = [-15 10 960 0 0 0 0 0 0]';
P0 = diag([10 10 10 10^4 10^4 10^4 5 5 5]);

[x_est, ~] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);

figure()
    hold on
    plot(x_est(1,:),'-', 'LineWidth', 2, 'DisplayName', 'p_x: X_{est}(1)')
    plot(x_est(2,:),'-', 'LineWidth', 2, 'DisplayName', 'p_y: X_{est}(2)')
    plot(x_est(3,:),'-', 'LineWidth', 2, 'DisplayName', 'p_z: X_{est}(3)')
    
    plot(p(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = p_x')
    plot(p(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = p_y')
    plot(p(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_3 = p_z')
    grid minor
    xlabel("Nro. iteraciones"); ylabel("estimaciones"); legend('show','FontSize', 15 )
    title("Estimacion con Kalman, x_{0/−1} = [−15 10 960 0 0 0 0 0 0]^T, P_{0/−1} = diag([10 10 10 10^4 10^4 10^4 5 5 5])")

%% (b) Para las posiciones, estados iniciales lejanos al real con varianza inicial peque˜na.
%x0/−1 = [10 − 600 50 0 0 0 0 0 0]T
%P0/−1 = diag([10 10 10 104 104 104 5 5 5])

X0 = [10 -600 50 0 0 0 0 0 0]';
P0 = diag([10 10 10 10^4 10^4 10^4 5 5 5]);

[x_est, ~] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);

figure()
    hold on
    plot(x_est(1,:),'-', 'LineWidth', 2, 'DisplayName', 'p_x: X_{est}(1)')
    plot(x_est(2,:),'-', 'LineWidth', 2, 'DisplayName', 'p_y: X_{est}(2)')
    plot(x_est(3,:),'-', 'LineWidth', 2, 'DisplayName', 'p_z: X_{est}(3)')
    
    plot(p(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = p_x')
    plot(p(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = p_y')
    plot(p(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_3 = p_z')
    grid minor
    xlabel("Nro. iteraciones"); ylabel("estimaciones"); legend('show','FontSize', 15 )
    title("Estimacion con Kalman, x_{0/−1} = [10 − 600 50 0 0 0 0 0 0]^T, P_{0/−1} = diag([10 10 10 10^4 10^4 10^4 5 5 5])")

%% (c) Para las posiciones, estados iniciales lejanos al real con varianza inicial alta.
%x0/−1 = [10 − 600 50 0 0 0 0 0 0]T
%P0/−1 = diag([106 106 106 104 104 104 5 5 5])

X0 = [10 -600 50 0 0 0 0 0 0]';
P0 = diag([10^6 10^6 10^6 10^4 10^4 10^4 5 5 5]);

[x_est, ~] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);

figure()
    hold on
    plot(x_est(1,:),'-', 'LineWidth', 2, 'DisplayName', 'p_x: X_{est}(1)')
    plot(x_est(2,:),'-', 'LineWidth', 2, 'DisplayName', 'p_y: X_{est}(2)')
    plot(x_est(3,:),'-', 'LineWidth', 2, 'DisplayName', 'p_z: X_{est}(3)')
    
    plot(p(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = p_x')
    plot(p(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = p_y')
    plot(p(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_3 = p_z')
    grid minor
    xlabel("Nro. iteraciones"); ylabel("estimaciones"); legend('show','FontSize', 15 )
    title("Estimacion con Kalman, x0/−1 = [10 − 600 50 0 0 0 0 0 0]^T, P0/−1 = diag([10^6 10^6 10^6 10^4 10^4 10^4 5 5 5])")
