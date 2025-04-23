% Generar una se�al de ejemplo
f=1000;
fs=2*f;
t = 0:1/fs:5; % Vector de tiempo [seg]
%signal = sin(2*pi*10*t); % Se�al sinusoidal de 10 Hz
rng(1);%para generar la misma se�al aleatoria, es para version de MATLAB R2014
signal = randn(length(t),1);

% Calcular la autocorrelaci�n del ruido "blanco"
[autocorr, lags] = xcorr(signal, 'unbiased');

% Especificaciones del filtro
orden = 4; % Orden del filtro
%frecuencia_corte = 100; % Frecuencia de corte en Hz
frecuencia_corte = 10; % Frecuencia de corte en Hz
frecuencia_muestreo = fs; % Frecuencia de muestreo en Hz

% Dise�o del filtro Butterworth
[numerador, denominador] = butter(orden, frecuencia_corte/(frecuencia_muestreo/2), 'low');

% Aplicaci�n del filtro a una se�al
X = filter(numerador, denominador, signal);%Se�al X
f=3; %[Hz]
Y = 0.5+0.5*square(2*pi*f*t,50);

% Calcular la autocorrelaci�n X
[autocorr_X, lags_X] = xcorr(X, 'unbiased');
% Calcular la autocorrelaci�n Y
[autocorr_Y, lags_Y] = xcorr(Y,'unbiased');color='r';

% Longitud de los datos
N = length(Y);

% Inicializar el vector de autocorrelaci�n
autocorr_Y_ = zeros(1, N);

% Calcular la media de los datos
mean_data = mean(Y);

% Calcular la autocorrelaci�n
for lag = 0:N-1
    sum_product = 0;
    for i = 1:N-lag
        sum_product = sum_product + (Y(i) - mean_data) * (Y(i + lag) - mean_data);
    end
    autocorr_Y_(lag + 1) = sum_product / (N - lag);
end

% Normalizar la autocorrelaci�n
autocorr_Y_ = autocorr_Y_ / autocorr_Y_(1);
N_=length(autocorr_Y_);
%autocorr_Y=[flip(autocorr_Y_) autocorr_Y_(1:N-1)];color='b';
%lags_Y=lag*2+1;

% Graficar la se�al original y su autocorrelaci�n
figure(1);
subplot(2,1,1);
plot(t, signal);
title('Ruido Blanco');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags, autocorr);
title('Autocorrelaci�n ruido Blanco');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelaci�n');

figure(2);
subplot(2,1,1);
plot(t, X);
title('X');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_X, autocorr_X);
title('Autocorrelaci�n X');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelaci�n X');

figure(3);
subplot(2,1,1);
plot(t, Y);
title('Y');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_Y,autocorr_Y,color);
title('Autocorrelaci�n Y');grid on;hold on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelaci�n Y');