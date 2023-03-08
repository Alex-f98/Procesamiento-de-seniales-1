clc; clear; close all;
%Problema 1
load('tp5.mat'); 
[N, M] = size(p);
%     a: [351×3 double]
%     p: [351×3 double]
%     v: [351×3 double]

%p: posicion; v: velocidad;  a: aceleracion
T = 1;      %[s]
%a'(t) = 0, [a'_x, a'_y, a'_z] ~ N(0, q)
% x = x0 + v.t + 0,5.a0.t^2 --> x = x0 + v.t / v = ∫a / a ~ N(0,q)
%a' = 0
%v' = a
%x' = v   -> X' = [x' v' a']^T = [0 1 0; 0 0 1; 0 0 0][x v a]' + [0 0 q]'
%Pero x, v, a son en realidad en 3D.
%x = [px py pz vx vy vz ax ay az]'
%Matriz de transicion de estados
Fc = [0 0 0, 1 0 0, 0 0 0;
     0 0 0, 0 1 0, 0 0 0;
     0 0 0, 0 0 1, 0 0 0;%     
     0 0 0, 0 0 0, 1 0 0;
     0 0 0, 0 0 0, 0 1 0;
     0 0 0, 0 0 0, 0 0 1;%
     0 0 0, 0 0 0, 0 0 0;
     0 0 0, 0 0 0, 0 0 0;
     0 0 0, 0 0 0, 0 0 0 ];
 
%q¯_c = [0,...,0,N(0,q)]T : ruido del proceso
q = diag([10^-2, 8e-3, 2e-5]);
      
%Qc = E{q¯c q¯c^T} : matriz de covarianza del ruido
Qc = zeros(size(Fc)); Qc(end-2:end, end-2:end) = q;


%% Problema 2
%a)
% F (nxn)(9,9)  G (nxm)(9,9)   H (lxn)(9,9)    R (lxl)(9,9)    Q (mxm)(9,9)  P (nxn)(9,9)   K (nxl)(9,9)
Fd = eye(size(Fc)) + Fc*T + (Fc*T)^2/2 ;      Gd = eye(size(Fc));
%Fd = eye(size(Fc)) + Fc*T + (Fc*T)^2/2 + (Fc*T)^3/3;      Gd = eye(size(Fc));

Qd = Qc*T;
%Qd = Qc*T + (Fc*Qc + 0.5*Qc*Fc')*T^2;

Hd = eye(9,9); %[mxn]

Sigma_p2 = 2500;    %[m^2]
Sigma_v2 = 25;      %[m^2/s^2]
Sigma_a2 = 1;       %[m^2/s^4]

%Rd = eye(M,M)*Sigma_p2;
Rd = diag([ones(1,3)*Sigma_p2, ones(1,3)*Sigma_v2, ones(1,3)*Sigma_a2]); %[mxm]

%-%b)
X0 = [10 -600 50 0 0 0 0 0 0]';
P0 = diag([10^6 10^6 10^6 10^4 10^4 10^4 5 5 5]);

%nu_k = sqrt(Sigma_p2)*randn(M, N);
nu_k = [sqrt(Sigma_p2)*randn(M, N);
        sqrt(Sigma_v2)*randn(M, N);
        sqrt(Sigma_a2)*randn(M, N)];


yk = Hd*[p, v, a]' + nu_k;

[x_est, innovaciones] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);
%x_est = [px py pz vx vy vz ax ay az]'

%%
figure()
    hold on
    plot(x_est(1,:),'LineWidth', 2, 'DisplayName', 'p_x: X_{est}(1)')
    plot(x_est(2,:),'LineWidth', 2, 'DisplayName', 'p_y: X_{est}(2)')
    plot(x_est(3,:),'LineWidth', 2, 'DisplayName', 'p_z: X_{est}(3)')
    
    plot(p(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = p_x')
    plot(p(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = p_y')
    plot(p(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_3 = p_z')
    grid minor
    xlabel("Nro. iteraciones"); ylabel("estimaciones"); legend('show','FontSize', 15 )
    title("Estimacion con Kalman")
 %%   
figure()
    hold on
    plot(x_est(4,:),'LineWidth', 2, 'DisplayName', 'v_x: X_{est}(4)')
    plot(x_est(5,:),'LineWidth', 2, 'DisplayName', 'v_y: X_{est}(5)')
    plot(x_est(6,:),'LineWidth', 2, 'DisplayName', 'v_z: X_{est}(6)')
    
    plot(v(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_4 = v_x')
    plot(v(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_5 = v_y')
    plot(v(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_6 = v_z')
    grid minor
    xlabel("Nro. iteraciones"); ylabel("estimaciones"); legend('show','FontSize', 15 )
    title("Estimacion con Kalman")
    
figure()
    hold on
    plot(x_est(7,:),'LineWidth', 2, 'DisplayName', 'v_x: X_{est}(7)')
    plot(x_est(8,:),'LineWidth', 2, 'DisplayName', 'v_y: X_{est}(8)')
    plot(x_est(9,:),'LineWidth', 2, 'DisplayName', 'v_z: X_{est}(9)')
    
    plot(a(:, 1), 'LineWidth', 2, 'DisplayName', 'Valor real: x_7 = a_x')
    plot(a(:, 2), 'LineWidth', 2, 'DisplayName', 'Valor real: x_8 = a_y')
    plot(a(:, 3), 'LineWidth', 2, 'DisplayName', 'Valor real: x_9 = a_z') 
    grid minor
    
    xlabel("Nro. iteraciones"); ylabel("estimaciones"); legend('show','FontSize', 15 )
    title("Estimacion con Kalman")
%%
figure()
%Nota: para el caso de medición de posición, incluya tambiín estas mediciones en el
%gráfico de la trayectoria para compararlas con las trayectorias real y estimada.
hold on
plot3(p(:,1), p(:,2), p(:,3))
plot3(x_est(1,:), x_est(2,:), x_est(3,:))
xlabel("p_x"); ylabel("p_y"); zlabel("p_z")
grid on
view(3)

%-%c) Calcule en cada uno de los casos anteriores la autocorrelaci´on de las innovaciones y verifique
%la validez del algoritmo de Kalman observando si dichas innovaciones son un proceso blanco.

figure()
    hold on;
    rx = xcorr(innovaciones(1,:));
    ry = xcorr(innovaciones(2,:));
    rz = xcorr(innovaciones(3,:));
    plot(rx, 'DisplayName', 'p_x'); plot(ry, 'DisplayName', 'p_y'); plot(rz, 'DisplayName', 'p_z')
    title("Autocorrelacion de las innovaciones , p")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")
    grid minor
    legend('show','FontSize', 15 )
    
figure()
    hold on;
    rx = xcorr(innovaciones(4,:));
    ry = xcorr(innovaciones(5,:));
    rz = xcorr(innovaciones(6,:));
    plot(rx, 'DisplayName', 'v_x'); plot(ry, 'DisplayName', 'v_y'); plot(rz, 'DisplayName', 'v_z')
    title("Autocorrelacion de las innovaciones , v", 'LineWidth', 3)
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")
    grid minor
    legend('show','FontSize', 15 )
    
figure()
    hold on;
    rx = xcorr(innovaciones(7,:));
    ry = xcorr(innovaciones(8,:));
    rz = xcorr(innovaciones(9,:));
    plot(rx, 'DisplayName', 'a_x'); plot(ry, 'DisplayName', 'a_y'); plot(rz, 'DisplayName', 'a_z')
    title("Autocorrelacion de las innovaciones , a")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")
    grid minor
    legend('show','FontSize', 15 )
    
%%
%d) Determine la observabilidad del sistema para los tipos de medición. Analice los resultados
% relacionándolos con la cantidad de estados que considere se estimaron correctamente.
Ob = obsv(Fd, Hd);
unobsv = length(Fd) - rank(Ob);
if unobsv; disp("No es observable");else; disp("Es observable"); end

