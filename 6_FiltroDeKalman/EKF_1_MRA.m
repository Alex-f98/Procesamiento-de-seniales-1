clc; clear; close all;
load('MRA.mat')

%k_e = 100;  %N/m
%M = 10;     %Kg
%b = 10;     %Kg/s
global M;
global b_;
global T;

syms 'x1' 'x2' 'x3' 'real'

%(error en el modelo con 𝜎b^2 = 0.2)
b_ = 18;     %Kg/s
sigma_b2 = 0.2;
T = 0.01;

%% a) asumiento que no hay fuerzas externas(u(t)=0): m.a(t) = -k_e.d(t) -b.v(t), hay cond. inicial
%x1 = v    x2=d (velocidad y posicion)

%Se quieren estimar los estados velocidad (x1) y posición (x2) 
%mediante kalman, pero suponiendo que se desconoce el valor de la constante elástica del resorte. Reformule el
%sistema de estados para poder estimar los estados

%    m.x2'= -x3.x1 - b.x2;
%f-> m.x2'+ x3.x1 + b.x2 = 0;

%f = @(x1,x2,x3)(m*x2'+ x3*x1 + b_*x2);

%x1,k+1 = (1 - T*b/m)*x1 - T*x2*x3/m + eps_b
%x2,k+1 = x1,k*T + x2,k
%x3,k+1 = x3,k

fun_F =@(x1,x2,x3) ([(1 - T*b/M)*x1 - T*x2*x3/M ; x1*T + x2; x3]);

%Fk = jacobian(f);
fun_Fk =@(x1,x2,x3) ([ (1 -b/M)*T, -x3*T/M, -T*x2/3;
                  T,        1,        0;
                  0,        0,        1]);
      
Hk = [0 1 0];

%Kalman extendido es un parche, no se demostro convergencia.


%% b) yk = dk + wk,  T = 10ms, wk~N(0; 0.05) es el ruido de midicion
sigmaW2 = 0.025;
%sigma_b2 = 0.2;                 %A MAYOR Q, CONFIO MAS EN LAS MEDICIONES QUE EN EL MODELO, POR ESO ES MAS RUIDOSO.

Gd = eye(3,3);                  %G (nxm)

%Qd = E[Vk.Vk^H]:
Q = diag([sigma_b2 0 0]);          %Cov(eps_v * eps_v')
Qd = Q*T;                       % Q (mxm)

%Rd = sigmaW2;
wk = sqrt(sigmaW2)*randn(1, length(v));

%a)
yk = d.' + wk;

P0 = diag([3 3 100].^2);
X0 = zeros(3, 1);
%b)Rd= 1,2;

Rd = sigmaW2;


%%
N = length(yk);
Y = yk;
 [l, ~] = size(Y);
n = length(X0);

Xk = X0;               %x_0/-1 = E[x0]
Pk = P0;               %P_0/-1 = Cov[x - x0
Rk = Rd;    
Gk = Gd;
Qk = Qd;
for k = 1:N
        
        hk_pri = Hk*Xk;
        
        %Actualizacion:
        K = Pk*Hk'/(Hk*Pk*Hk' + Rk);     %Kk = P_k/k-1 H'_k/(H_k*P_k/k-1*H'_k + R_k)
        %ek = Y(:,k) - Hk*Xk;            %ek = yk - H.X_k/k-1
        ek = Y(:,k) - hk_pri;             %ek = yk - H.X_k/k-1
        
        Xk = Xk + K*ek;               %Tiene que ser un proceso Blanco.
        Pk = (eye(n,n) - K*Hk)*Pk;    %P_k/k = (I - K_k*H_k)*P_k/k-1 

        %Almaceno datos:
        innovaciones(:, k) = ek;
        x_est(:, k) = Xk;
        Ks(:, k) = K;

        fk_post = fun_F(Xk(1), Xk(2), Xk(3));
        Fk = fun_Fk(Xk(1), Xk(2), Xk(3));
        
        %Prediccion:
        %Xk = Fk*Xk;                 %X_k+1/k = F_k*X_k/k
        Xk = fk_post;
        Pk = Fk*Pk*Fk' + Gk*Qk*Gk'; %P_k+1/k = F_k*P_k/k*F'_k + G_k.Q_k.G'_k


end

%%


figure()
    hold on
    plot(x_est(1,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(1)')
    plot(x_est(2,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(2)')
    plot(x_est(3,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada: X_{est}(3)')
    
    plot(yk, 'DisplayName', 'señal observada: y_k')
    
    plot(d, 'LineWidth', 2, 'DisplayName', 'Valor real: x_2 = d(n)')
    plot(v, 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = v(n)')
    plot([1, length(v)],[Ke Ke], 'LineWidth', 2, 'DisplayName', 'Valor real: x_3 = Ke')
    
    grid minor
    title("Estimacion con Kalman, b = 18 \neq 10, Q = 0")
    xlabel("Nro. iteraciones")
    ylabel("estimaciones")
    legend show



