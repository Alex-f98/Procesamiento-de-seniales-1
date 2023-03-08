clc;
clear;
close all;

%Filtro IIR
%Se requiere diseñar un filtro discreto IIR, pasa bajos, tipo Butterworth.
%Especificaciones
delta_p = 0.2;
delta_s = 0.1;
wp = 0.32*pi;
ws = 0.6*pi;
T = 2;

%Transformación bilineal
Omega_p = (2/T)*tan(wp/2);
Omega_s = (2/T)*tan(ws/2);

d = (((1-delta_p)^(-2) - 1)/(delta_s^(-2) - 1))^(1/2);
N = log10(d)/log10(Omega_p/Omega_s);

N = ceil(N);

Omega_c_max = Omega_s/((delta_s)^(-2) - 1)^(1/(2*N));
Omega_c_min = Omega_p/((1-delta_p)^(-2) - 1)^(1/(2*N));

Omega_c = (Omega_c_min + Omega_c_max)/2;

k = 0:1:2*N-1;

%Polos Butterworth Pasa Bajos
s = Omega_c*exp(1i*((N+1+2*k)*pi/(2*N)));
s = s(1:N);

%Discretización del filtro diseñado
z = (1 + (T/2)*s)./(1 - (T/2)*s);

A_k = s./(s - T/2);
A = prod(A_k);

s_k = s;
z_k = z;

%Transferencia en el dominio continuo
num = prod(-s_k);
den = poly(s_k);
H = tf(num,den);
pzmap(H)
grid on

%Transferencia en el dominio discreto 

%Términos del numerador de la transferencia
%Los ceros de H(Z) son N -1s.
b = A*poly(-1*ones(N,1));

%Términos del denominador de la transferencia
%Los polos de H(z) son los zk.
a = poly(z_k);

[H,w] = freqz(b,a,w);
figure
plot(w/pi,abs(H),'DisplayName','Transferencia del filtro diseñado')
grid on;
grid minor;
xlabel('$\frac{\omega}{\pi}$ [rad/s]','interpreter','latex','FontSize',16)
ylabel('|Hk(w)|','FontSize',16)
legend()

figure
zplane(b,a)
grid
title('Diagrama de polos y ceros')

