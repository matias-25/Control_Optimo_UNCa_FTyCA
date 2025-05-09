clear all, close all;
%% Parametros
f = 3;
T = 1/f;
fs = 100*f;
% *******Tiempo final en funcion del periodo T********
t = 0:1/fs:5*T; % Vector de tiempo [seg]
N = length(t);
A = 2;
sigma = 0.6;
wo = 2*pi*f;
%señal x
x = A.*sin(wo*t-sigma);

%% Calcular la autocorrelación x
[autocorr_x, lags_x] = xcorr(x,x, 'unbiased');

%% Calcular la autocorrelación x sin la funcion xcorr
integral =0;
tau =0;
x_tau = zeros(1,N);
phi_xx_se = zeros(1,N);
t_tau = 0:1/fs:5*T;

%t_tau=[-flip(t(2:N)) t];
k=1;
figure(1);subplot(2,1,1);plot(t,x,'r');title(' x  vs x_\tau');xlabel('tiempo [seg]');
hold on;grid on;
tau_selec=0.5*(N-1)*1/fs;% cual tau se quiere plotear
for tau=0:1/fs:5*T
    x_tau = A.*sin(wo*(t-tau)-sigma);
    if tau==tau_selec
    subplot(2,1,1);plot(t_tau,x_tau);hold on;legend('x','x_\tau');
    %tau_selec=0;
    end
        for j=1:N
            integral = integral + x(j)*x_tau(j);
        end
    phi_xx_se(k)= integral;
    
    k=k+1;

end
phi_xx = (max(x)/2)*((phi_xx_se/max(phi_xx_se))+1-0.0322);
%phi_xx =phi_xx_se;
subplot(2,1,2);
plot(t_tau,phi_xx,'c');title('Autocorrelación x');xlabel('\tau [seg]');

%% Graficar la señal X y su autocorrelación
figure(2);
subplot(2,1,1);
plot(t, x);
title('x');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_x(N:2*N-1)/(N-1), autocorr_x(N:2*N-1));
title('Autocorrelación x');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación x');