%%
function [x_est, innovaciones] = Kalman( Y, P0, Qk, X0, Fk, Gk ,Hk, Rk)

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
    %Ks = zeros(n, N);
   
    
    for k = 1:N

        %Actualizacion:
        K = Pk*Hk'/(Hk*Pk*Hk' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/k-1*H'_k + R_k)
        %K = Pk*Hk'*pinv(Hk*Pk*Hk' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/k-1*H'_k + R_k)
        ek = Y(:,k) - Hk*Xk;            %ek = yk - H.X_k/k-1)
        Xk = Xk + K*ek;               %Tiene que ser un proceso Blanco.
        Pk = (eye(n,n) - K*Hk)*Pk;    %P_k/k = (I - K_k*H_k)*P_k/k-1 

        %Almaceno datos:
        innovaciones(:, k) = ek;
        x_est(:, k) = Xk;
        %Ks(:, k) = K;


        %Prediccion:
        Xk = Fk*Xk;                 %X_k+1/k = F_k*X_k/k
        Pk = Fk*Pk*Fk' + Gk*Qk*Gk'; %P_k+1/k = F_k*P_k/k*F'_k + G_k.Q_k.G'_k

        %if k == N-1;disp("P_{end}:" ); Pk; end
    end
end
