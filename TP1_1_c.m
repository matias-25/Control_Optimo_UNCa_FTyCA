clear all;
%sistema a lazo abierto
A=[-2  2;2 -5];
B=[1;0];
C=[0 1];
D=0;
%sistema ampliado
A_amp=[A zeros(2,1);-C 0];
B_amp=[B;0];

%Matriz Q:Pondera los estados del sistema
%Q=diag([1 1 10]);color='b';
Q=diag([1 10 1]);color='c';
%Q=diag([1 1 1]);color='r';
%Matriz R:Pondera la señal de control
%R=0.1;
R=1;
K_amp=lqr(A_amp,B_amp,Q,R);
polos=eig(A_amp-B_amp*K_amp)
tao_min=abs(min(polos));%tiempo de muestreo
tao_max=abs(max(polos));
At=1/(10*tao_min);
T_max=5*1/tao_max;
N_max=round(T_max/At);
X=[0;0];
ref=1;
psi=0;
y=zeros(1,N_max);
accion=y;
e_amp=[X;psi];
for k=1:N_max
	   
    u=-K_amp*[X;psi];
    xp=A*X+B*u;
    psi=psi+At*(ref-C*X);
    X=X+At*xp;
    y(k)=C*X;
    accion(k)=u;
end
tiempo=0:At:T_max;
figure(1);
subplot(2,1,1);
plot(tiempo(1:length(y)),y,color);xlabel('tiempo [seg]');title('Salida h2');grid on;hold on;
subplot(2,1,2);
plot(tiempo(1:length(y)),accion,color);xlabel('tiempo [seg]');title('Accion de Control');grid on;hold on;