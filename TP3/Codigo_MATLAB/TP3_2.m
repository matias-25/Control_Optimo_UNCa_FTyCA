clear all, close all;
%% Parametros
f = 10;
fs = 100*f;
t = 0:1/fs:1; % Vector de tiempo [seg]
N = length(t);
A = 2;
sigma = 1;
wo = 2*pi*f;
%señal x
x = A.*sin(wo*t-sigma);

%% Calcular la autocorrelación x
[autocorr_x, lags_x] = xcorr(x, 'unbiased');

%% Graficar la señal X y su autocorrelación
figure(1);
subplot(2,1,1);
plot(t, x);
title('x');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_x(N:2*N-1), autocorr_x(N:2*N-1));
title('Autocorrelación x');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación x');