clear all;
%sistema a lazo abierto
A=[-2  2;2 -5];
B=[1;0];
C=[0 1];
D=0;
figure(1);
sys=ss(A,B,C,D);
rlocus(sys);

%obtension de la ecuacion caraterística
%|sI-A|=[s 0;  -  [-2  2;
    %    0 s]       2 -5]
%|sI-A|= |s+2 -2; -2 s+5|=(s+2)(s+5)-(-2)*(-2)=s^2+5s+2s-4+10=s^2+7s+6  
%p=poly(A);
%matriz de polos deseados
J=[-6+1i -6-1i];

%Matriz de controlabilidad
%M=ctrb(A,B);
%det_M=det(M);

%Diseño mediante la Fórmula de Ackermann
Kack=acker(A,B,J);
K=place(A,B,J)

% Calcular la matriz del sistema en lazo cerrado
Acl = A - B * K;

% Calcular los polos del sistema en lazo cerrado
poles = eig(Acl)
h=1/(10*max(abs(poles)));
%Diseño control óptimo en tiempo continuo con un integrador en el lazo para 
%referencia distinta de cero
tiempo=10/h;
t=0:h:tiempo*h;
X=[0 ; 0];
A_amp=[A zeros(2,1);-C 0];
B_amp=[B;0];
C_amp=[0 1 0];
J_amp=[-6+1i -6-1i -1];color='b';
%J_amp=[-6+1i -6-1i -10];color='r';
K_amp=place(A_amp,B_amp,J_amp)
i=0;
ref=1; %%altura en metros del tanque 2
psi=0; %Error integrado

while(i<(tiempo))
  
  u=-K_amp*([X;psi]);%ley de control
  X_P=A*X+B*u;%X punto
  X=X+h*X_P;%Esto es el cálculo de la integral como sumatoria
  psi=psi+h*(ref-C*X);
  
  i=i+1;
  h1(i)=X(1);
  h2(i)=X(2);
  accion(i)=u;
end

figure(2);
subplot(2,1,1);
plot(t,h2,color); title('Salida del sistema');xlabel('tiempo[seg]');grid on; hold on;
%step(sys,'r'), legend('altura controlada','altura sin controlador');
subplot(2,1,2) ;
plot(t,accion,color); title('Accion de control "caudal"');xlabel('tiempo[seg]');grid on;hold on;
