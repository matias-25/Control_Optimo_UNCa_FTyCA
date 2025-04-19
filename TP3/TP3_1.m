% Generar una señal de ejemplo
t = 0:0.001:1; % Vector de tiempo [seg]
%signal = sin(2*pi*10*t); % Señal sinusoidal de 10 Hz
rng(1);%para generar la misma señal aleatoria, es para version de MATLAB R2014
signal = randn(length(t),1);

% Calcular la autocorrelación
[autocorr, lags] = xcorr(signal, 'coeff');

% Graficar la señal original y su autocorrelación
figure;
subplot(2,1,1);
plot(t, signal);
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags, autocorr);
title('Autocorrelación');
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación');
