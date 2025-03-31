clear all;close all;
%sistema a lazo abierto
A=[-2  2;
    2 -5];
B=[1;0];
C=[0 1];
D=0;
[num , den]= ss2tf(A,B,C,D);
sys= tf(num,den);
f=step(sys);
h=0.01;
f1 = diff(f);%derivada primera
f2 = diff(f1);%derivada segunda
p_inflex_t=0;
p_inflex=0;
pend_tg=0;
t=0:h:(length(f)*h)-h;

%punto de inflexion
for i=1:length(f)-2
    if i<length(f)-2;
    if (f2(i)>0 && f2(i+1)<0) || (f2(i)<0 && f2(i+1)>0)
        p_inflex_t=t(i);
        p_inflex=f(i);
        pend_tg=(f(i)-f(i-1))/(t(i)-t(i-1));
    end
    end
end
rec_tg=pend_tg*(t-p_inflex_t)+p_inflex;%recta tangente al punto de inflexion

figure(1);
plot(t,f);title('Planta a Lazo Abierto');xlabel('tiempo[seg]');axis([0 5 -0.2 1.5]);grid on; hold on;
plot(p_inflex_t,p_inflex,'O');hold on
plot(t,rec_tg,'r');

%obtencion de L y T_ para calcular los parámetros del controlador PID
%con el primer metodo de Ziegler-Nichols
for i=1:length(f)
    if 0.995<rec_tg(i) && 1.005>rec_tg(i)
        T_=t(i);
    end
    if min(abs(rec_tg(i)))<=0.005 
       L=t(i); 
    end
end
%Las constantes PID
Kp=1.2*T_/L;
Ti=2*L;Ki=Kp/Ti;
Td=0.5*L;Kd=Kp*Td;
%PID
A1=((2*Kp*h)+(Ki*(h^2))+(2*Kd))/(2*h);
B1=(-2*Kp*h+Ki*(h^2)-4*Kd)/(2*h);
C1=Kd/h;
e=zeros(length(f),1);%error
u=0;
X=-[0;0];
i=0;k=0;
referencia=1; %%altura en metros del tanque 2
for j=0:h:(length(f)*h)-h;
  i=i+1;
  k=i+2;
  e(k)=referencia - X(2); %ERROR
  u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
  accion(i)=u;
  h2(i)=X(2);
  for jj=0:h
    X_P=A*X+B*u;
    X=X+h*X_P;
  end
end

figure(2);
subplot(2,1,1);
plot(t,accion);title('Accion');xlabel('tiempo[seg]');grid on;
subplot(2,1,2);
plot(t,h2);title('Salida del sistema');xlabel('tiempo[seg]');grid on;
