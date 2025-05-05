clear all, close all;
f=1000;
fs=2*f;
t = 0:1/fs:1; % Vector de tiempo [seg]
rng(1);%para generar la misma señal aleatoria, es para version de MATLAB R2014
signal = randn(length(t),1);

%% Calcular la autocorrelación del ruido "blanco"
[autocorr, lags] = xcorr(signal, 'coeff');

%% Especificaciones del filtro
orden = 4; % Orden del filtro
%frecuencia_corte = 100; % Frecuencia de corte en Hz
frecuencia_corte = 10; % Frecuencia de corte en Hz
frecuencia_muestreo = fs; % Frecuencia de muestreo en Hz
%% Diseño del filtro Butterworth
[numerador, denominador] = butter(orden, frecuencia_corte/(frecuencia_muestreo/2), 'low');

%% Señales X,Y,Z,W
X = filter(numerador, denominador, signal);%Señal X
f_Y=3; %[Hz]
Y = 0.5+0.5*square(2*pi*f_Y*t,50);
Z =1 * ones(size(t));
W = zeros(1,length(t));
    for i=1:length(t)
        W(i)= X(i)+Y(i)+Z(i);
    end

%% Calcular la autocorrelación X
[autocorr_X, lags_X] = xcorr(X, 'coeff');

%% Calcular la autocorrelación Y
[autocorr_Y, lags_Y] = xcorr(Y,'unbiased');color='r';

%% Calcular la autocorrelación Z
[autocorr_Z, lags_Z] = xcorr(Z,'unbiased');

%% Calcular la autocorrelación W
[autocorr_W, lags_W] = xcorr(W,'unbiased');

%% Graficar la señal original y su autocorrelación
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

%% Graficar la señal X y su autocorrelación
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

%% Graficar la señal Y y su autocorrelación
figure(3);
subplot(2,1,1);
plot(t, Y);
title('Y');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_Y,autocorr_Y,color);
title('Autocorrelación Y');grid on;hold on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación Y');

%% Graficar la señal Z y su autocorrelación
figure(4);
subplot(2,1,1);
plot(t, Z);
title('Z');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_Z,autocorr_Z);
title('Autocorrelación Z');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación Z');

%% Graficar la señal W y su autocorrelación
figure(5);
subplot(2,1,1);
plot(t, W);
title('Z');grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(lags_W,autocorr_W);
title('Autocorrelación W');grid on;
xlabel('Lags');
ylabel('Coeficiente de Autocorrelación W');