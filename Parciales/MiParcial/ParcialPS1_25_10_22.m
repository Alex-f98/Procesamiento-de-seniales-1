clc; clear; close all;
load('x.mat'); % se carga en una variable 'x'
N = length(x);
L = N/2; %12; %N/4;
%n = 0:N;
%K debe ser #frec x 2
%se sabe que el menor error posible esta cerca de la mitad de N.
%en el osea L ~ N/2
xn = x;
%% a)
yy = 0;
for j = 1: N-L+1
   yL = transpose(flip(xn(j:j+L-1)));
   yy = yy + yL*yL';
end
R_ = yy/N;                      %matriz de correlacion LxL

% Calcular la descomposición en autovalores y autovectores
% R = [S G]Γ [S^H; G^H]
% donde S contiene los primeros K autovectores.
[V, LAMBDA] = eig(R_);
[V_ordenado, idx]= sort(diag(LAMBDA), 'descend');

figure()
plot(V_ordenado, '*r')
xlim([0 50])
grid on
%En x = 7 ya se ve como los autovalores son aproximadamente los mismos esto
%quiere decir que desde 7 en adelante pertenecen  a S.
%entonces K = 6.
K = 6;

%% b)
w = SPRIT( L, K, xn);
%Como K = 6 y las frecuencias dada por el enunciado son:
%   0,45*pi .. 0.54*pi .. 0.1035*pi (3 frecuencias) entonces las otras
%   faltantes seran negativas.
w_ordenado = sort(w/pi, 'descend');
%tomo los primeros 3 de mayor a manor.
w_ref = [0.54 0.45 0.1035];

Errores = [];
for L = K+1: N-5
    w = SPRIT( L, K, xn);
    w_ordenado = sort(w/pi, 'descend');
    ErrorPorc = 100* abs(w_ordenado(1:K/2) - w_ref')./(w_ref');
    Errores(end+1) = max(ErrorPorc);
end

figure()
hold on
plot([0 N-5], [1 1], '-.k')
plot(K+1:1:N-5, Errores, '-b', 'LineWidth', 2)
xlabel("L")
xlim([20 100])
ylabel("Error relativo porcentual [%]")
grid on

%Aca se ve que desdde L = 53 el error es menor al 1%, luego hasta caso el
%final donde todo aumenta abruptamente.
%%
%verifiquemos con L = 12
%L = ceil(N/2);
L = 53;
w = SPRIT( L, K, xn);
w_ordenado = sort(w/pi, 'descend');
%tomo los primeros 3 de mayor a manor.
w_ref = [0.54 0.45 0.1035];


w_ordenado = sort(w/pi, 'descend');
ErrorPorc = 100* abs(w_ordenado(1:K/2) - w_ref')./(w_ref');
disp("W"+w_ordenado+"pi")
disp("Error porcentual w: " + ErrorPorc + "%")

%% sprit
function [w] = SPRIT( L, K, xn)
%a)
N = length(xn);
%L = 15; %estimar K.
%K = 6;  %segun las frecuencias de yn x2
% A partir de y[1], · · · y[N], estimar la matriz de covarianza R_.
yy = 0;
for j = 1: N-L+1
   yL = transpose(flip(xn(j:j+L-1)));
   yy = yy + yL*yL';
end
R_ = yy/N;                      %matriz de correlacion LxL

% Calcular la descomposición en autovalores y autovectores
% R = [S G]Γ [S^H; G^H]
% donde S contiene los primeros K autovectores.
[V, LAMBDA] = eig(R_);
[~, idx]= sort(diag(LAMBDA), 'descend');
V = V(:,idx);
S = V(:, 1:K);                               %[LxK]

% Obtener S1 y S2 quitando la última y la primera fila de S
S1 = S(1:L-1, :);           %[0 eye(L-1)]*S;
S2 = S(2:L,   :);           %[0 eye(L-1)]*S;

% Calcular los autovalores de [S^H_1 S_1]^−1  [S^H_1 S2_]
[~, Ds] = eig( inv(S1'*S1) * (S1'*S2) );

% Obtener las frecuencias ωi a partir de los autovalores anteriores.
w = angle(diag(Ds));
%disp("frecuencia/pi:"+ w/pi );
end
