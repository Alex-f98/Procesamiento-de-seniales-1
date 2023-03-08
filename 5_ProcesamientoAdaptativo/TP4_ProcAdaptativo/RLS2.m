%RLS2
function [Weigths, errors, e_V, x_estimado] = RLS2(u, x, v, M, delta, lambda)
%[x, u, v, var_v] = fun(sn, SNR_dB);
    L = length(u);
    errors = zeros(1,L);
    x_estimado = zeros(1,L);
    e_V = zeros(1, L);
    %Inicializo el algoritmo RLS
    Weigths = zeros(M,L-M+1);
    Weigths(:,1) = [3;4];
    Pn = (1/delta)*eye(M);
    %Kn = zeros(M,1);

    for j = 1:L-M+1
        yn = flip(u(j:j+M-1));
        Kn = lambda\Pn*yn/(1 + lambda\yn'*Pn*yn);       %Ganancia(Mx1)
                                                        x_est = Weigths(:,j)'*yn;
        error = x(j+M-1) - x_est;                       %Error a priori 
                                                        x_estimado(j) = x_est;
                                                        errors(j+M-1) = error;
                                                        e_V(j+M-1) = (v(j+M-1)-x_est);
        Weigths(:,j+1) = Weigths(:,j) + Kn*error';      %Actualizacion
        Pn = (1/lambda)*Pn - (1/lambda)*Kn*yn'*Pn;      %Mat inv de correlacion
    end
   
    
end
    