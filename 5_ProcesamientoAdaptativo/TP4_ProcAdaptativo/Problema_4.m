clc; clear; close all
% En este problema se buscan comparar algunas características entre los algoritmos adaptativos LMS y RLS. 
% En particular se analizarán las diferencias de convergencia y mínimo error
% entre ambos. Considere en todos los casos M = 2 y w(0) = [3 4 ]T
M = 2;
w_inicial = [3, 4];
SNR_dB = 10;
L = 4000;
sigma_g2 = 6.42e-4;

%% a) lambda = 0.99, delta = 0.001, m>=200 J(n)_ V(n)_ y sigma_g2
delta = 0.001;
lambda = 0.99;
mu_LMS = 40; 

%Número de realizaciones
m = 500;

VV_ = zeros(m,L);
JJ_ = zeros(m,L);

for i=1:m

    gn = sqrt(sigma_g2)*randn(L, 1);
    h = [1 0.5 0.1 0.3 0.4 0.24];
    sn = filter(h, 1, gn);

    [xn, un, vn, var_v] = fun(sn, SNR_dB);
    
    [w_rls, errors, e_V, x_est] = RLS2(un, xn, vn, M, delta, lambda);
    JJ_(i,:) = errors.^2;
    VV_(i,:) = e_V.^2;
    
end

%Curva de aprendizaje
J = mean(VV_, 1);

%Curva de la diferencia al cuadrado entre el ruido y la salida del filtro
V = mean(JJ_, 1);

%figuras:-------------
figure()
hold on
semilogy(1:L, J, '-r', 'LineWidth', 2)%, 'DisplayName', "J")
semilogy(1:L, V, '-b', 'LineWidth', 2)%, 'DisplayName', "V")
semilogy([1 L], [var(sn) var(sn)], '-k')

title("Curva de aprendizaje con 500 iteraciones, \lambda = 0.99, \delta = 0.001")
xlabel("Nro de iteraciones")
ylabel("$\hat{J}(n)$, $\hat{V}(n)$",'Interpreter','latex')
grid minor

lgd = legend('$\hat{J}(n)$','$\hat{V}(n)$','$\sigma_{s}^2$');
set(lgd,'Interpreter','latex'),set(lgd,'FontSize',12);
set(lgd,'Location','northeast');

%% c) 
% Implemente ambos algoritmos, RLS y LMS, considerando λ = 0,998 y modificando el paso
% µ de LMS para que ambos tengan aproximadamente la misma velocidad de convergencia
% (buscando la misma pendiente inicial de las potencias de error). Grafique V(n)_ para ambos
% casos y justifique las diferencias observadas.


delta = 0.001;
lambda = 0.998;
mu_LMS = 470; 

%Número de realizaciones
m = 500;

VV_RLS = zeros(m,L);    JJ_RLS_ = zeros(m,L);
VV_LMS = zeros(m,L);    JJ_LMS_ = zeros(m,L);

for i=1:m

    gn = sqrt(sigma_g2)*randn(L, 1);
    h = [1 0.5 0.1 0.3 0.4 0.24];
    sn = filter(h, 1, gn);

    [xn, un, vn, var_v] = fun(sn, SNR_dB);
    
    [Weigths_RLS, errors_RLS, e_V_RLS, x_est_RLS] = RLS2(un, xn, vn, M, delta, lambda);
    [Weigths_LMS, errors_LMS, e_V_LMS, x_est_LMS] = LMS2(un, xn, vn,  M, mu_LMS, w_inicial);
    
    JJ_RLS(i,:) = errors_RLS.^2;
    VV_RLS(i,:) = e_V_RLS.^2;
    
    JJ_LMS(i,:) = errors_LMS.^2;
    VV_LMS(i,:) = e_V_LMS.^2;
    
end

%Curva de aprendizaje
J_RLS = mean(VV_RLS, 1);

%Curva de la diferencia al cuadrado entre el ruido y la salida del filtro
V_RLS = mean(JJ_RLS, 1);

%Curva de aprendizaje
J_LMS = mean(VV_LMS, 1);

%Curva de la diferencia al cuadrado entre el ruido y la salida del filtro
V_LMS = mean(JJ_LMS, 1);


%figuras:-------------
figure()
hold on
semilogy(1:L, J_RLS, '-r', 'LineWidth', 2)%, 'DisplayName', "J")
semilogy(1:L, V_RLS, '-b', 'LineWidth', 2)%, 'DisplayName', "V")

semilogy(1:L, J_LMS, '-g', 'LineWidth', 2)%, 'DisplayName', "J")
semilogy(1:L, V_LMS, '-m', 'LineWidth', 2)%, 'DisplayName', "V")

semilogy([1 L], [var(sn) var(sn)], '-k')

title("Curva de aprendizaje con 500 iteraciones, \lambda = 0.998, \delta = 0.001, \mu = 470 ")
xlabel("Nro de iteraciones")
ylabel("$\hat{J}(n)$, $\hat{V}(n)$",'Interpreter','latex')
grid minor

lgd = legend('$\hat{J}(n)_{RLS}$','$\hat{V}(n)_{RLS}$', ...
    '$\hat{J}(n)_{LMS}$','$\hat{V}(n)_{LMS}$','$\sigma_{s}^2$');
set(lgd,'Interpreter','latex'),set(lgd,'FontSize',17);
set(lgd,'Location','northeast');

%% c) 
% Nuevamente corra ambos algoritmos, manteniendo λ = 0,998, pero esta vez ajustando
% el par´ametro µ para que la potencia Vb(n) de ambos tienda al mismo valor m´ınimo para
% n → ∞. Grafique ambas curvas y analice las diferencias observadas.

delta = 0.001;
lambda = 0.998;
mu_LMS = 27; 

%Número de realizaciones
m = 500;

VV_ = zeros(m,L);
JJ_ = zeros(m,L);

for i=1:m

    gn = sqrt(sigma_g2)*randn(L, 1);
    h = [1 0.5 0.1 0.3 0.4 0.24];
    sn = filter(h, 1, gn);

    [xn, un, vn, var_v] = fun(sn, SNR_dB);
    
    [Weigths_RLS, errors_RLS, e_V_RLS, x_est_RLS] = RLS2(un, xn, vn, M, delta, lambda);
    [Weigths_LMS, errors_LMS, e_V_LMS, x_est_LMS] = LMS2(un, xn, vn,  M, mu_LMS, w_inicial);
    
    JJ_RLS(i,:) = errors_RLS.^2;
    VV_RLS(i,:) = e_V_RLS.^2;
    
    JJ_LMS(i,:) = errors_LMS.^2;
    VV_LMS(i,:) = e_V_LMS.^2;
    
end

%Curva de aprendizaje
J_RLS = mean(VV_RLS, 1);

%Curva de la diferencia al cuadrado entre el ruido y la salida del filtro
V_RLS = mean(JJ_RLS, 1);

%Curva de aprendizaje
J_LMS = mean(VV_LMS, 1);

%Curva de la diferencia al cuadrado entre el ruido y la salida del filtro
V_LMS = mean(JJ_LMS, 1);


%figuras:-------------
figure()
hold on
semilogy(1:L, J_RLS, '-r', 'LineWidth', 2)%, 'DisplayName', "J")
semilogy(1:L, V_RLS, '-b', 'LineWidth', 2)%, 'DisplayName', "V")

semilogy(1:L, J_LMS, '-g', 'LineWidth', 2)%, 'DisplayName', "J")
semilogy(1:L, V_LMS, '-m', 'LineWidth', 2)%, 'DisplayName', "V")

semilogy([1 L], [var(sn) var(sn)], '-k')

title("Curva de aprendizaje con 500 iteraciones, \lambda = 0.998, \delta = 0.001, \mu = 27 ")
xlabel("Nro de iteraciones")
ylabel("$\hat{J}(n)$, $\hat{V}(n)$",'Interpreter','latex')
grid minor

lgd = legend('$\hat{J}(n)_{RLS}$','$\hat{V}(n)_{RLS}$', ...
    '$\hat{J}(n)_{LMS}$','$\hat{V}(n)_{LMS}$','$\sigma_{s}^2$');
set(lgd,'Interpreter','latex'),set(lgd,'FontSize',17);
set(lgd,'Location','northeast');

