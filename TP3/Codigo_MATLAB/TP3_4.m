clear all; close all;
%%  Densidad espectral de potencia
K=1;T=10;
Fw=tf([0 K],[T 1]) %Transformada de Fourier

%Grafico de Bode
figure(1); 
bode(Fw);title('Densidad espectral de potencia S(\omega)'); grid on;

%%  Función de correlación 
t=0:0.01:5; %vector tiempo
% la transformada inversa de la Densidad espectral de potencia S(w) 
phi_xx = K/T*exp(-t/T);

% Analizar la influencia de T en la función de correlación
% se incrementará y decrementará T exponencialmente
M = 5;
figure(2); 
%incremento exponencial de T
subplot(2,1,1);
title('Funcion de correlacion, K=1, incremento exponencial de T');
color= 'r';hold on; grid on;
xlabel('tiempo [s]'); ylabel('amplitud');
for i=1:M
   switch i
    case 1
           color = 'r';
    case 2
           color = 'm';
    case 3
           color = 'c';
    case 4
           color = 'b';
    case 5
           color = 'g';
   end
   T = 2^(i);
   plot(t,K/T*exp(-t/T), color);hold on; grid on;
end
legend('T=2^1','T=2^2','T=2^3','T=2^4','T=2^5');

%decremento exponencial de T
subplot(2,1,2);
title('Funcion de correlacion, K=1, decremento exponencial de T');
color= 'b';hold on; grid on;
xlabel('tiempo [s]'); ylabel('amplitud');
for i=1:M
   switch i
    case 1
           color = 'r';
    case 2
           color = 'm';
    case 3
           color = 'c';
    case 4
           color = 'b';
    case 5
           color = 'g';
   end
   T = 2^(-i);
   plot(t,K/T*exp(-t/T), color);hold on; grid on;  
end
legend('T=2^-^1','T=2^-^2','T=2^-^3','T=2^-^4','T=2^-^5');
