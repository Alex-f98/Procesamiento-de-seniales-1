%%
function [x_est, innovaciones, Ks] = EKF( Y, P0, Qk, X0, fk, Gk ,hk, Rk)

    syms 'x1' 'x2' 'x3' 'real'
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

    Fk_ = @(x1, x2, x3)(jacobian(fk, [x1, x2, x3]));
    Hk_ = @(x1, x2, x3)(jacobian(hk,[x1, x2, x3]));
    
        Hk = feval(Hk_, [X0(1),X0(2),X0(3)]');
        Fk = feval(Fk_, [X0(1),X0(2),X0(3)]');
    
    [~, n] = size(Fk);  [~, m] = size(Gk);   [l, ~] = size(Hk);

    Xk = X0;               %x_0/-1 = E[x0]
    Pk = P0;               %P_0/-1 = Cov[x - x0]

    innovaciones = zeros(l, N);
    x_est = zeros(n, N);
    %Ks = zeros(n, N);
   
    
    for k = 1:N
        x1 = Xk(1);
        x2 = Xk(2);
        x3 = Xk(3);
        
        Hk = feval(Hk_, [x1,x2,x3]');
        Fk = feval(Fk_, [x1,x2,x3]');
        
        %Actualizacion:
        K = Pk*Hk'/(Hk*Pk*Hk' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/k-1*H'_k + R_k)
        %K = Pk*Hk'*pinv(Hk*Pk*Hk' + Rk);  %Kk = P_k/k-1 H'_k/(H_k*P_k/k-1*H'_k + R_k)
        %ek = Y(:,k) - Hk*Xk;            %ek = yk - H.X_k/k-1)
        ek = Y(:,k) - feval(hk, [x1,x2,x3]');            %ek = yk - H.X_k/k-1)
        
        Xk = Xk + K*ek;               %Tiene que ser un proceso Blanco.
        Pk = (eye(n,n) - K*Hk)*Pk;    %P_k/k = (I - K_k*H_k)*P_k/k-1 

        %Almaceno datos:
        innovaciones(:, k) = ek;
        x_est(:, k) = Xk;
        Ks(:, k) = K;


        %Prediccion:
        %Xk = Fk*Xk;                 %X_k+1/k = F_k*X_k/k
        Xk = feval(fk, [x1,x2,x3]');
        Pk = Fk*Pk*Fk' + Gk*Qk*Gk'; %P_k+1/k = F_k*P_k/k*F'_k + G_k.Q_k.G'_k


    end
end
