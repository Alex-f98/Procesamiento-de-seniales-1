clc; clear; close all
% En este problema se busca analizar el comportamiento del algoritmo LMS en función de
% diferentes parámetros. Para poder visualizar mejor los resultados con un mayor número de
% realizaciones independientes, se deberá utilizar una señal sintética en lugar del audio. Para
% ello, defina el código para generar la señal como un proceso s(n) de largo L = 4000 cuya
% entrada g(n) sea ruido blanco gaussiano de media nula y varianza Sigmag^2 = 6,42 × 10−4
% tal que
% 
% s(n) = g(n) + 0,5 g(n − 1) + 0,1 g(n − 2) + 0,3 g(n − 3) + 0,4 g(n − 4) + 0,24 g(n − 5).

L = 4000;
sigma_g2 = 6.42e-4;
gn = sqrt(sigma_g2)*randn(L, 1);

h = [1 0.5 0.1 0.3 0.4 0.24];
sn = filter(h, 1, gn);

%% a) Curva de aprendizaje:
mu_LMS = 40; 
M = 2;
w_inicial = [3; 4];
SNR_dB = 10;
MM = 500;

% [xn, un, vn, varR] = fun(sn, SNR_dB);
% [Weigths, errors, x_est] = LMS(un, xn, M, mu_LMS, w_inicial);


JJ_ = zeros(MM, L);
Vn_ = zeros(MM, L);

for i = 1:MM
    
    gn = sqrt(sigma_g2)*randn(L, 1);
    h = [1 0.5 0.1 0.3 0.4 0.24];
    sn = filter(h, 1, gn);

    [xn, un, vn, varR] = fun(sn, SNR_dB);

    %[Weigths, errors, x_est] = LMS(un, xn, M, mu_LMS, w_inicial);
    [Weigths_LMS, errors_LMS, e_V_LMS, x_est_LMS] = LMS2(un, xn, vn,  M, mu_LMS, w_inicial);
    
    JJ_(i,:) = errors_LMS.^2;
    Vn_(i,:) = e_V_LMS.^2 ;
    
end

% Ry = zeros(M,M);
% for i = M : L; Ry = Ry + un(i:-1:i-M+1).*un(i:-1:i-M+1)'/L; end
% disp("2/Traza(R_y):"+ (2/trace(Ry)));

%Curva de aprendizaje.
J = mean(JJ_, 1); % la media de cada columna.
V = mean(Vn_, 1);

%figuras:-------------
figure()
hold on
plot(1:L, J, '-r')%, 'DisplayName', "J")
plot(1:L, V, '-b')%, 'DisplayName', "V")
plot([1 L], [var(sn) var(sn)], '-k')

title("Curva de aprendizaje con 500 iteraciones")
xlabel("Nro de iteraciones")
ylabel("$\hat{J}(n)$, $\hat{V}(n)$",'Interpreter','latex')
grid minor

lgd = legend('$\hat{J}(n)$','$\hat{V}(n)$','$\sigma_{s}^2$');
set(lgd,'Interpreter','latex'),set(lgd,'FontSize',12);
set(lgd,'Location','northeast');

%% b) Parametrización con M:
mu_LMS = 50;

for M = [2 3 10 30]
     
    %JJ_ = zeros(MM, L-1);
    w_inicial = zeros(M,1);
    Vn_ = zeros(MM, L);

    for i = 1:MM

        gn = sqrt(sigma_g2)*randn(L, 1);
        h = [1 0.5 0.1 0.3 0.4 0.24];
        sn = filter(h, 1, gn);

        [xn, un, vn, varR] = fun(sn, SNR_dB);

        %[Weigths, errors, x_est] = LMS(un, xn, M, mu_LMS, w_inicial);
        [~, ~, e_V_LMS, ~] = LMS2(un, xn, vn,  M, mu_LMS, w_inicial);
    
        Vn_(i,:) = e_V_LMS.^2;
        %JJ_(i,:) =  ;

    end

    %Curva de aprendizaje.
    %J = mean(JJ_, 1); % la media de cada columna.
    V = mean(Vn_, 1);
    
    %figuras:-------------
    figure(5)
    hold on
    semilogy(1:L, V, 'DisplayName', sprintf('M = %g',M)) 
    %plot(1:L-1, J, 'r', 'DisplayName', sprintf('M = %g',M)) 
end
    legend show
     
    title("Curva de aprendizaje con 500 iteraciones")
        xlabel("Nro de iteraciones")
        ylabel("$\hat{V}(n)$",'Interpreter','latex')
        grid minor

%% c) Parametrización con µ:

M = 2;
for mu_LMS = [50, 160, 500, 1600]
     
    %JJ_ = zeros(MM, L-1);
    w_inicial = zeros(M,1);
    Vn_ = zeros(MM, L);

    for i = 1:MM

        gn = sqrt(sigma_g2)*randn(L, 1);
        h = [1 0.5 0.1 0.3 0.4 0.24];
        sn = filter(h, 1, gn);

        [xn, un, vn, varR] = fun(sn, SNR_dB);

       % [Weigths, errors, x_est] = LMS(un, xn, M, mu_LMS, w_inicial);
        [~, errors_LMS, e_V_LMS, ~] = LMS2(un, xn, vn,  M, mu_LMS, w_inicial);
    
        %JJ_(i,:) = errors_LMS.^2;
        Vn_(i,:) = e_V_LMS.^2 ;

    end

    %Curva de aprendizaje.
    %J = mean(JJ_, 1); % la media de cada columna.
    V = mean(Vn_, 1);
    
    %figuras:-------------
    figure(10)
    hold on
    plot(1:L, V, 'DisplayName', sprintf('\\mu = %g',mu_LMS)) 
    %plot(1:L-1, J, 'r', 'DisplayName', sprintf('M = %g',M)) 
end

    legend show
     
    title("Curva de aprendizaje con 500 iteraciones")
    xlabel("Nro de iteraciones")
    ylabel("$\hat{V}(n)$",'Interpreter','latex')
    grid minor
