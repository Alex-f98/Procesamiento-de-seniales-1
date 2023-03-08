function [Weigths, errors, x_est] = LMS(yn, x, N, u, w_inicial) 
%yn(vector columna)
%w_inicial (vector columna)
assert( (size(yn,2)==1 || size(yn,2)==1),"ERR: x(n) e y(n) tienen que ser vectores columna");

if ~exist('w_inicial')
        w_inicial = zeros(N,1); 
end
if size(w_inicial,2)~=1; w_inicial = w_inicial'; end
Weigths =[];  
errors = [];
x_est =[];

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
   for i = 1:L-N+1
        Y_    = flip(yn(i:i+N-1));
        xn_   = w' * Y_;
        error = x(i+N-1) - xn_;
        w     = w + u * Y_ * error; %conj(error);
        
        Weigths(:, end+1) = w;
        errors(end+1) = error;
        x_est(end+1) = xn_; 
   end
   
end