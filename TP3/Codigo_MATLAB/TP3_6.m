%% señal secuencia numérica pseudoaleatoria de 1000 datos en el intervalo [0,1)
N=1000000; %numeros de datos
x=0;xv=zeros(N,1);

for hh=1:N
    x=((pi+x)^5)-fix((pi+x)^5);
    xv(hh)=x;
end

%% Autocorrelacion de la señal
%hist(xv)
%rng(1);
%xv=rand(N,1);
fixx = xcorr(xv, 'unbiased'); %no sesgada
fixx_nor = xcorr(xv, 'coeff'); %normalizada


figure(1);
subplot(2,1,1);
plot(fixx(N:2*N-1));title('Autocorrelación de x');
xlabel('Tiempo \tau');
subplot(2,1,2);
plot(fixx_nor(N:2*N-1));title('Autocorrelación de x Normalizado');
xlabel('Tiempo \tau');

%% Espectro de potencia Sxx
SxM=fft(fixx(N:2*N-1));
SxM_nor=fft(fixx_nor(N:2*N-1));
figure(2);
subplot(2,1,1);
semilogx(20*log10(abs(SxM(1:N/2))),'.k');title('Módulo de F(w) en dB'); grid on;
xlabel('Pulsación en rad por seg');
ylabel('S_{xx}(w)');
subplot(2,1,2);
semilogx(20*log10(abs(SxM_nor(1:N/2))),'.k');title('Módulo de F(w) en dB con autocorrelacion normalizada'); grid on;
xlabel('Pulsación en rad por seg');
ylabel('S_{xx}(w)');
