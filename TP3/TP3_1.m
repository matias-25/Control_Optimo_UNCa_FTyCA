% Generar una señal de ejemplo
f=1000;
fs=2*f;
t = 0:1/fs:pi; % Vector de tiempo [seg]
%signal = sin(2*pi*10*t); % Señal sinusoidal de 10 Hz
rng(1);%para generar la misma señal aleatoria, es para version de MATLAB R2014
signal = randn(length(t),1);

% Calcular la autocorrelación del ruido "blanco"
[autocorr, lags] = xcorr(signal, 'coeff');

% Especificaciones del filtro
orden = 4; % Orden del filtro
%frecuencia_corte = 100; % Frecuencia de corte en Hz
frecuencia_corte = 10; % Frecuencia de corte en Hz
frecuencia_muestreo = fs; % Frecuencia de muestreo en Hz

% Diseño del filtro Butterworth
[numerador, denominador] = butter(orden, frecuencia_corte/(frecuencia_muestreo/2), 'low');

% Aplicación del filtro a una señal
X = filter(numerador, denominador, signal);%Señal X
Y = 0.5+0.5*square(3*2*pi*t,50);
%Y = square(2*pi*t,50);
% Calcular la autocorrelación X
[autocorr_X, lags_X] = xcorr(X, 'coeff');
% Calcular la autocorrelación Y
[autocorr_Y, lags_Y] = xcorr(Y, 'coeff');
% Graficar la señal original y su autocorrelación
figure(1);
subplot(2,1,1);
plot(t, signal);
title('Ruido Blanco');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags, autocorr);
title('Autocorrelación ruido Blanco');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación');

figure(2);
subplot(2,1,1);
plot(t, X);
title('X');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_X, autocorr_X);
title('Autocorrelación X');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación X');

figure;
subplot(2,1,1);
plot(t, Y);
title('Y');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_Y, autocorr_Y);
title('Autocorrelación Y');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación Y');