clear all; close all;
%% Generación de ruido blanco
fs = 1000; t_max = 10; sigma = 1;
t = 0:1/fs:t_max; N = length(t);

x = sigma * randn(1, N);
figure(1);
subplot(3, 1, 1);
plot(t,x);title('Ruido Blanco');xlabel('t[seg]'); ylabel('Amplitud');

%% Autocorrelacion 
[phi_xx, lags] = xcorr(x, 'unbiased'); tau = lags / fs;
subplot(3, 1, 2);
plot(tau, phi_xx, 'b-'); xlabel('Desplazamiento \tau (s)'); ylabel('\phi_{xx}(\tau)'); title('Función de Autocorrelación del Ruido Blanco'); grid on; xlim([0 0.5]);


%% Autoespectro de Potencia
S_xx = fft(phi_xx); f = (-N/2:N/2-1)*(fs/N);
subplot(3, 1, 3);
semilogx(20*log10(abs(S_xx(1:N/2))),'.k');
xlabel('Frecuencia (Hz)'); ylabel('S_{xx}(f)'); title('Autoespectro de Potencia del Ruido Blanco'); grid on; xlim([-fs/2 fs/2]);