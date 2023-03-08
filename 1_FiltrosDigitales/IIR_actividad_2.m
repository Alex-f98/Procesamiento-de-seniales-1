clc;
clear;
close all;

%Se requiere diseñar un filtro discreto, pasa altos, IIR tipo Butterworth, con las
%especificaciones que se indican en la figura.

nfft = 16384;

%Especificaciones
ws = 0.2*pi;
wp = 0.42*pi;
delta_s = 0.1;
delta_p = 0.22;
T = 1;

%Transformación bilineal
%Especificaciones para un filtro pasa altos
Omega_p_pa = (2/T)*tan(wp/2);
Omega_s_pa = (2/T)*tan(ws/2);

%Especificaciones para un filtro pasa bajos
Omega_p_pb = 1/Omega_p_pa;
Omega_s_pb = 1/Omega_s_pa;

%Cálculo del orden del filtro
d = (((1-delta_p)^(-2) - 1)/(delta_s^(-2) - 1))^(1/2);
N = log10(d)/log10(Omega_p_pb/Omega_s_pb);
N = ceil(N)

%Frecuencia de corte aproximada para el filtro pasa bajos
Omega_c_max = Omega_s_pb/((delta_s)^(-2) - 1)^(1/(2*N));
Omega_c_min = Omega_p_pb/((1-delta_p)^(-2) - 1)^(1/(2*N));
Omega_c = (Omega_c_min + Omega_c_max)/2

%Polos Butterworth Pasa Bajos
k = 0:1:2*N-1;
s = Omega_c*exp(1i*((N+1+2*k)*pi/(2*N)));
s_k = s(1:N);

%Polos Butterworth Pasa Altos
s_k_pa = 1./s_k;

%Discretización del filtro diseñado
z_k_pa = (1 + (T/2)*s_k_pa)./(1 - (T/2)*s_k_pa);

A_k = 1./(1 - (T/2)*s_k_pa);
A = prod(A_k)

%Transferencia en el dominio continuo del filtro PA
num = [1 zeros(1,N-1)];
den = poly(s_k_pa);
H = tf(num,den);
pzmap(H)
grid on

%Transferencia en el dominio discreto del filtro PA

%Términos del numerador de la transferencia
%Los ceros de H(Z) son N 1s.
b = A*poly(1*ones(N,1))

%Términos del denominador de la transferencia
%Los polos de H(z) son los zk.
a = poly(z_k_pa)

%[h,w] = freqz(___,n,'whole') devuelve la respuesta de frecuencia en n puntos 
%de muestra alrededor de todo el círculo de la unidad

[H,w] = freqz(b,a,nfft,'whole');
figure
plot(w/pi,abs(H),'DisplayName','Transferencia del filtro diseñado')
grid on;
grid minor;
xlabel('$\frac{\omega}{\pi}$ [rad/s]','interpreter','latex','FontSize',16)
ylabel('|H(w)|','FontSize',16)
legend()

figure
zplane(b,a)
grid
title('Diagrama de polos y ceros del filtro PA')

