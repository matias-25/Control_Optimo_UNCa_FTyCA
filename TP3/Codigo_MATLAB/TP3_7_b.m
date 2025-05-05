clear all;close all;
N=500;m=4;
y=zeros(N,1);
x=zeros(1,m+1);x(1)=1;
  for k=1:N
    res=xor(x(4),x(5)); %PRBS4 x^4+x^3+1
    y(k)=2*res-1;
    x_d=circshift(x,[1,1]);
    x_d(1)=res;
    x=x_d;
  end
figure
subplot(3,1,1);
plot(y);title('PRBS4');
fixx = xcorr(y,y, 'unbiased');
subplot(3,1,2);
plot(fixx(N:2*N-1));title('Autocorrelación de x');
xlabel('Tiempo tao')
SxM=fft(fixx(N:2*N-1));
subplot(3,1,3);
semilogx(20*log10(abs(SxM(1:N/2))),'.k');
xlabel('Pulsación en rad por seg')
title('Módulo de F(w) en dB')