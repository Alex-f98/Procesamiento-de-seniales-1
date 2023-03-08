function [Weigths, errors, e_V, x_est] = LMS2(yn, x, v, M, u, w_inicial) 
%yn(vector columna)
%w_inicial (vector columna)
L = length(yn);
assert( (size(yn,2)==1 || size(yn,2)==1),"ERR: x(n) e y(n) tienen que ser vectores columna");

if ~exist('w_inicial')
        w_inicial = zeros(M,1); 
end
if size(w_inicial,2)~=1; w_inicial = w_inicial'; end
% Weigths =[];  
% errors = [];
% x_est =[];

Weigths = zeros(M,L-M+1);
Weigths(:,1) = w_inicial;

errors = zeros(1,L);
x_est= zeros(1,L);
e_V = zeros(1, L);
%u; %tiene un rango 0 < mu < 1/traza(A)
%a) ALGORITMO LMS
L = length(yn);
% ➢ definir las condiciones iniciales W(0)
    w = w_inicial;
% ➢ Paso 0: estimar W(1) a partir de W(0)
% ➢ Paso n: estimar W(n+1) a partir de W(n):
    % 1. Se calcula la salida del filtro,          x_(n)= W'_n.Y(n).
    % 2. Se calcula el error de estimación,        e(n) = x(n) - x_(n).
    % 3. Se adaptan los coeficientes del filtro,   W{n+1}_= Wn_ + mu.Y(n)e(n)*
   for i = 1:L-M+1
        Y_    = flip(yn(i:i+M-1));
        xn_   = w' * Y_;
        error = x(i+M-1) - xn_;
        w     = w + u * Y_ * error; %conj(error);
        Weigths(:,i+1) = Weigths(:,i) + u * Y_ * error;
        
                                                        x_est(i) = xn_;
                                                        errors(i+M-1) = error;
                                                        e_V(i+M-1) = (v(i+M-1)-xn_);
        
%         Weigths(:, end+1) = w;
%         errors(end+1) = error;
%         x_est(end+1) = xn_; 
   end
   
end