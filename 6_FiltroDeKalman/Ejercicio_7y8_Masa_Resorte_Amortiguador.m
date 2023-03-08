%Una secuencia de datos binarios se transmiten a traves de un canal de
%comunicaciones con respuesta impusiva hk y ruido blanco gaussiano de media
%nula y varianza sigma_W2 a la salida, tal que 
% yk = h0.s_k + h1.s_k-1 + h2.s_k-2, como se indica en la fig.
%Suponiendo que se conocen los coeficientes del canal, se desea
%implementar, se un ecualizador que permite recuperar la secuencia de datos
%s_k mediante un filtro de Kalman. Para ello se dispone de las
%observaciones z_k = yk + wk a la salida del sistema.
%                      wk--.
%                           \
%       sk ----| H(z) |---->(+)--->zk

%(a) definir F y H el ruido del proceso de donde viene?
%z_k = yk + wk       yk = h0.s_k + h1.s_k-1 + h2.s_k-2
%z_k = h0.s_k + h1.s_k-1 + h2.s_k-2  +  wk
%
%  ---Señal-->z_k = [h0 h1 h2][s_k s_k-1 s_k-2]' + wk
%  x1,k = s_k ,       x2,k = s_k-1  ,    x3,k = s_k-2     (estados)
%  x1,k+1 = s_k+1 ,   x2,k+1 = x1,k  ,   x3,k+1 = x2,k    (estados derivada)
%                       X_deriv = [x1,k+1 x2,k+1 x3,k+1]
%
%        X(k+1) =       F               X(k) +      Vk 
%   [X_deriv]^T = [0 0 0; 1 0 0; 0 1 0][X]^T + [s_k+1 0 0]^T.
%El ruido del proceso proviene de una entrada del sistema.
Fd = [0 0 0;
      1 0 0;
      0 1 0];
  
Hd = [0.4 0.3 0.2];


%b) sk = 2.bk -1 / bk ~ Bernoulli(0.5), N=1000, hk = [0,4 0,3 0,2],
%deltaW2=10^-4
N = 1000;
bk = binornd(1, 0.5, [1 N]);
sk = 2*bk - 1;
sigma_w2 = 10^-4;
sigma_s2 = var(sk);

wk = sqrt(sigma_w2)*randn(1,N);

X0 = [0; 0; 0];
P0 = diag([10 10 10]);

X = [sk; [0, sk(1:end-1)]; [0 0, sk(1:end-2)]];   %[s_k s_k-1 s_k-2]'
yk = Hd*X + wk;

Gd = eye(3);
Qd = [sigma_s2 0 0;
             0 0 0;
             0 0 0 ];
         
Rd = sigma_w2;

[x_est, innovaciones] = Kalman( yk, P0, Qd, X0, Fd, Gd ,Hd, Rd);


figure()
    hold on
    plot(x_est(1,:),'-o', 'LineWidth', 1, 'DisplayName', 'Señal estimada FK: X_{est}(1)')

    plot(yk, '-*', 'DisplayName', 'señal observada: y_k=H.X+w_k')
   
    plot(sk, 'LineWidth', 2, 'DisplayName', 'Valor real: x_1 = s_k')
   
    grid minor
    title("Estimacion con Kalman ")
    xlabel("Nro. iteraciones")
    ylabel("estimaciones")
    legend show
  
    
figure()
    r = xcorr(innovaciones);
    plot(r)
    title("Autocorrelacion de las innovaciones")
    xlabel("Nro. iteraciones")
    ylabel("Autorrelacion")




