% funcion RLS

function [w_, errors, Weigths, x_est] = RLS(y, xn, M, delta, lambda)
L = length(y);
errors = zeros(1,L-M+1);
x_est = zeros(1,L-M+1);
Weigths = zeros(M, L-M+1);
wn_ = zeros(M,1);                                        %w0_  0
Pn = eye(M)/delta;                                       %P0 = delta^-1.I_MxM
%Resto de las iteraciones: n: 1, 2,..., L-1.
for n = 1:L-M
   yn = flip(y(n:n+M-1));
   
   Kn = lambda\Pn*yn/(1 + lambda\yn'*Pn*yn);             %Ganancia(Mx1)
   eps = xn(n+M-1) - wn_'*yn;                            %Error a priori  
   x_est(n) = wn_'*yn;
   wn_ = wn_ + Kn * conj(eps);                           %Actualizacion
   Pn = lambda\Pn - lambda\Kn*yn'*Pn;                    %Mat inv de correlacion
   
   
   errors(n) = eps;
   Weigths(:, n+1) = wn_; %el primer peso estará asignado.
   
   
end
    w_ = wn_; 
end
