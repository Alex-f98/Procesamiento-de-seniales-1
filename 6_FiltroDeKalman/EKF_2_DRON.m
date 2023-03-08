clc; clear; close all;
global p1x p1y p1z p2x p2y p2z p3x p3y p3z p4x p4y p4z p5x p5y p5z p6x p6y p6z p7x p7y p7z p8x p8y p8z
%Problema 1
load('kalman_loc.mat')
[N, M] = size(p);
%     a: [351×3 double]
%     p: [351×3 double]
%     v: [351×3 double]

%p: posicion; v: velocidad;  a: aceleracion
T = 1;      %[s]
C = 3e8;

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

%Fd = eye(size(Fc)) + Fc*T + (Fc*T)^2/2 ;      Gd = eye(size(Fc));
Fd = eye(size(Fc)) + Fc*T + (Fc*T)^2/2 + (Fc*T)^3/3;      Gd = eye(size(Fc));

%Qd = Qc*T;
Qd = Qc*T + (Fc*Qc + 0.5*Qc*Fc')*T^2;

%asi para cada piy, es un lio, mejor definir el f y sacar de ahi

Sigma_p2 = 2500;    %[m^2]
Sigma_v2 = 25;      %[m^2/s^2]
Sigma_a2 = 1;       %[m^2/s^4]

Rd = eye(8,8)*Sigma_p2;

nu_k = sqrt(varw)*randn(8, N);
% 
% p1x = rp(1,1); p1y = rp(1,2); p1z = rp(1,3);
% p2x = rp(2,1); p2y = rp(2,2); p2z = rp(2,3);
% p3x = rp(3,1); p3y = rp(3,2); p3z = rp(3,3); 
% p4x = rp(4,1); p4y = rp(4,2); p4z = rp(4,3); 
% p5x = rp(5,1); p5y = rp(5,2); p5z = rp(5,3); 
% p6x = rp(6,1); p6y = rp(6,2); p6z = rp(6,3); 
% p7x = rp(7,1); p7y = rp(7,2); p7z = rp(7,3); 
% p8x = rp(8,1); p8y = rp(8,2); p8z = rp(8,3);
% 
% 
% yk_ = @(px, py, pz)([   sqrt( (px-p1x)^2 + (py-p1y)^2 + (pz - p1z)^2  )/C;
%                         sqrt( (px-p2x)^2 + (py-p2y)^2 + (pz - p2z)^2  )/C;
%                         sqrt( (px-p3x)^2 + (py-p3y)^2 + (pz - p3z)^2  )/C;
%                         sqrt( (px-p4x)^2 + (py-p4y)^2 + (pz - p4z)^2  )/C;
%                         sqrt( (px-p5x)^2 + (py-p5y)^2 + (pz - p5z)^2  )/C;
%                         sqrt( (px-p6x)^2 + (py-p6y)^2 + (pz - p6z)^2  )/C;
%                         sqrt( (px-p7x)^2 + (py-p7y)^2 + (pz - p7z)^2  )/C;
%                         sqrt( (px-p8x)^2 + (py-p8y)^2 + (pz - p8z)^2  )/C]);
% yk = zeros(8, N);
% for idx = 1:N
%     yk(:, idx) = yk_(p(idx,1), p(idx,2), p(idx,3));
% end
% yk = yk + nu_k;

yk = tau' + nu_k;

X0 = [10 ; -600 ; 50 ; 2 ; 3 ; 4 ; 0 ; 0 ; 0];
P0 = diag([10^5 10^5 10^5 10^3 10^3 10^3 0.9 0.9 0.9]);

%%

Y = yk;
%P0, 
Qk = Qd;
%X0, 
Fk = Fd;
Gk = Gd;
%Hk=  varia
Rk = Rd;

%function [x_est, innovaciones, Ks] = Kalman( Y, P0, Qk, X0, Fk, Gk ,Hk, Rk)

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
    n = length(X0);
    [l, ~] = size(Y);
    


    Xk = X0;               %x_0/-1 = E[x0]
    Pk = P0;               %P_0/-1 = Cov[x - x0]

    innovaciones = zeros(l, N);
    x_est = zeros(n, N);
    %Ks = zeros(n, N);
   
    h = zeros(8,1);
    Hk = zeros(8,9);
    for k = 1:N
        
        
        for idx = 1:length(rp)
            %h()
           h(idx) = ( (Xk(1)- rp(idx, 1))^2 + (Xk(2) - rp(idx, 2))^2 + (Xk(3) - rp(idx, 3))^2 )/C;     
           %Hk = d h(x)/ dx
           Hk(idx, :) = [ Xk(1) - rp(idx, 1), Xk(2) - rp(idx, 2), Xk(3) - rp(idx, 3), zeros(1,6)]./( h(idx)*C^2 );
        end
        

        %Actualizacion:
        K = Pk*Hk'/(Hk*Pk*Hk' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/k-1*H'_k + R_k)
        %K = Pk*Hk'*pinv(Hk*Pk*Hk' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/k-1*H'_k + R_k)
        ek = Y(:,k) - Hk*Xk;            %ek = yk - H.X_k/k-1)
        Xk = Xk + K*ek;                 %Tiene que ser un proceso Blanco.
        Pk = (eye(n,n) - K*Hk)*Pk;      %P_k/k = (I - K_k*H_k)*P_k/k-1 

        %Almaceno datos:
        innovaciones(:, k) = ek;
        x_est(:, k) = Xk;
        %Ks(:, k) = K;

        %Prediccion:
        Xk = Fk*Xk;                 %X_k+1/k = f(X_k/k)
        
        Pk = Fk*Pk*Fk' + Gk*Qk*Gk'; %P_k+1/k = F_k*P_k/k*F'_k + G_k.Q_k.G'_k


    end
%end

figure()
    hold on
    plot(x_est(1,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(1)')
    plot(x_est(2,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(2)')
    plot(x_est(3,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(3)')
    
    %plot(yk, 'DisplayName', 'señal observada: y_k')
    
    plot(p(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = p_x')
    plot(p(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = p_y')
    plot(p(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_3 = p_z')
    
    grid minor
    title("Estimacion con Kalman Extendido")
    xlabel("Nro. iteraciones")
    ylabel("estimaciones")
    legend show

