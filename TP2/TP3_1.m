% Generar una se�al de ejemplo
t = 0:0.001:1; % Vector de tiempo [seg]
%signal = sin(2*pi*10*t); % Se�al sinusoidal de 10 Hz
rng(1);%para generar la misma se�al aleatoria, es para version de MATLAB R2014
signal = randn(length(t),1);

% Calcular la autocorrelaci�n
[autocorr, lags] = xcorr(signal, 'coeff');

% Graficar la se�al original y su autocorrelaci�n
figure;
subplot(2,1,1);
plot(t, signal);
title('Se�al Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags, autocorr);
title('Autocorrelaci�n');
xlabel('Lags');
ylabel('Coeficiente de Autocorrelaci�n');