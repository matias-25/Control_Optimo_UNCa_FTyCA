clear all;close all;
%% Se define f(t) original
%% Parametros
f = 4; %frecuencia de la señal
fs = 1000*f; %frecuencia de muestreo
t_F = 1; %tiempo final en [seg]
t = 0:1/fs:t_F; % Vector de tiempo [seg]
A = 3; %Amplitud
T = 1/f; % Periodo
d = 0.5*T; %desplazamiento

%señal f
f_t = (A/2)+(A/2)*square(2*pi*f*(t+(d/2)),50);
%Grafico
figure(1);
plot(t,f_t); title('f(t) vs f reconstruida por Serie de Fourier, con primera armónica'); 
xlabel('tiempo en [seg]'); hold on;

%% Espectro de Fourier
t_1=0; t_2=T;
t_aux = t_1:1/fs:t_2;
w0 = 2*pi/T;
N = 5; %Numeros de armónicos
F = zeros(1,N+1);
vector_frecuencia = zeros(1,N+1); 
for n =0:N
    
    for i =1:length(t_aux)
     F(n+1) = F(n+1) + f_t(i)*exp(-(n)*w0*1i*t(i));
    end
    
  vector_frecuencia(n+1)= n;  
end
modudo_Fn= abs(F)/max(abs(F));
f_aux = 0:1/fs:N; % Crear el vector f_aux
Sinc_f = abs(sinc(d * f .* f_aux)); % Calcular Sinc_f
%Grafico Módulo de F
figure(2);
plot(vector_frecuencia, modudo_Fn ,'*r');hold on; title('Espectro en frecuencias  F_n'); xlabel('frecuencia');
plot(f_aux,Sinc_f, 'b');
%% Serie de Fourier
% Como exp(n*w0*1i*t)= cos(n*w0*t)+isen(n*w0*t)
% Entonces f(t) = F(n)*(cos(n*w0*t)+isen(n*w0*t));
% Suma de cos y sen reales : f(t) = (abs(F(n))/N)*(cos(n*w0*t)+sen(n*w0*t));
f_t_F = zeros(N,length(t)); % armónicos
f_t_rec = zeros(1,length(t)); % f(t) reconstruida
c = zeros(1,N);

for n =0:N
    c(n+1) = (A/2)*F(n+1)/max(abs(F));
    f_t_F(n+1,1:length(t)) = 2*c(n+1)*(cos((n)*w0*t)+1i*sin((n)*w0*t));
end
f_t_F(1,1:length(t)) = (f_t_F(1,1:length(t)))/2;
for n =0:N
f_t_rec = f_t_rec + f_t_F(n+1,1:length(t));
end

%Grafico
figure(3);
for j =1:N+1
    plot(t,f_t_F(j,1:length(t)));title('Armónicos '); xlabel('frecuencia');hold on;
end
figure(1);
plot(t,f_t_rec, 'r');plot(t,c(1)+ f_t_F(2,1:length(t)),'m');legend('f(t) original','f(t) Reconstruida','1° Armónico');
