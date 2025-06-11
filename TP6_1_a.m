clear all;%close all; 
%% Constantes del sistema
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5; 

%Condiciones iniciales
d(1)=0;
d_p(1)=0;
fi(1)=pi/3;
fi_p(1)=0;
u(1)=0; %accion de control 
x= [-d(1);-d_p(1);fi(1);fi_p(1)];
x0=x;

%Versión linealizada en el equilibrio inestable. Sontag Pp 104. 
Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0;0 0 0 1;0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_Bc=[0; 1/M; 0; -1/(long*M)];
Mat_Cc=[0 0 1 0; 1 0 0 0];

%% Sistema continuo
sys_c=ss(Mat_Ac,Mat_Bc,Mat_Cc,0); 

%% Sistema discreto
fs=100; %frecuencia de muestreo
Ts=1/fs;

sys_d=c2d(sys_c,Ts,'zoh'); 
Mat_Ad=sys_d.a;
Mat_Bd=sys_d.b;
Mat_Cd=sys_d.c;

%% controlador DLQG
Q=diag([1e-3 1e0 1e-4 1e0]);R=1e0;%Matrices de diseño del controlador DLQG
Kx = dlqr(Mat_Ad,Mat_Bd,Q,R); %ganacia del controlador

%% Obsevador de Luenberger
%Cálculo del Observador--------------------------------------------------- 
A_o=Mat_Ad'; 
B_o=Mat_Cd'; 
C_o=Mat_Bd'; 
Qo=diag([1e0 1e0 1e0 1e0]);Ro=diag([1e-3 1e0]); 
Ko= dlqr(A_o,B_o,Qo,Ro); %ganancia del observador

%% Monte Carlo
% Consigna 0, 0.01, 0.02, 0.05 y 0.1
sQ=0; %Para F
sR=0; %Para G. Covarianza del ruido de medicion sigma=sqrt(sR)
F_=sQ*eye(4); %Covarianza del ruido de estado Sigma=sqrt(sQ)
G_=sR;
S=Q;
P=S; %condición inicial de P
kmax=2000;
Realizaciones=5; %Cantidad de realizaciones para el Monte Carlo.
%Kx=zeros(kmax,4);
Kv=zeros(kmax,4);
Aa=Mat_Ad;
Ba=Mat_Bd;
u=zeros(Realizaciones,kmax);
Jn_=zeros(Realizaciones,kmax);
Jmin=x0'*P*x0;
J=0;
t=0:kmax-1;
y_sal=zeros(Realizaciones,kmax);

for hi=kmax-1:-1:1
    P= Q + Aa'*P*Aa - Aa'*P*Ba*inv(R+Ba'*P*Ba)*Ba'*P*Aa;
    %Kx(hi,:)=inv(R+Ba'*P*Ba)*Ba'*P*Aa;
    Kv(hi,:)=inv(R+Ba'*P*Ba)*Ba'*P*F_;
    %Ea(:,hi)=eig(Aa-Ba*Kx(hi,:));
end

for trial=1:Realizaciones %Empieza el Monte Carlo
    v=randn(4,kmax);%Señales aleatorios de media nula y varianza unidad.
    w=randn(2,kmax);
    x=x0+F_*v(:,1);
    x_hat=[0;0;0;0];% variables estado observador 
    d(trial,1)=x(1);
    d_p(trial,1)=x(2);
    fi(trial,1)=x(3);
    fi_p(trial,1)=x(4);
    
    for ki=1:kmax-1
        %u(trial,ki)=-Kx(1,:)*x; %DLQR
        %Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*Q*x + u(trial,ki)'*R*u(trial,ki));
        u(trial,ki)=-Kx(1,:)*x_hat-Kv(1,:)*v(:,ki); %LQG Ec 6-45,con observador
        Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*Qo*x + u(trial,ki)'*R*u(trial,ki));%Ec 11-38 con observador
        %y_sal(trial,ki)=Mat_Cd*x+G_*w(ki);
        Ys=Mat_Cd*x + G_*w(:,ki);
        y_sal(trial,ki)=Ys(1);
        x=mopdm2(Ts,x,u(trial,ki))+F_*v(:,ki+1);
        x_hat= Mat_Ad*x_hat + Mat_Bd*u(trial,ki) + Ko'*(Ys-Mat_Cd*x_hat);% Obsevador de Luenberger
        d(trial,ki+1)=x(1);
        d_p(trial,ki+1)=x(2);
        fi(trial,ki+1)=x(3);
        fi_p(trial,ki+1)=x(4);
    end
    Jn_(trial,ki+1)=Jn_(trial,ki+1)+x'*S*x;
    u(trial,ki+1)=-Kx(1,:)*[x]-Kv(1,:)*[v(:,ki)];
end

%% Graficos
t=t*Ts;
TamanioFuente=14;
figure;
subplot(3,2,2);hold on;grid on;title('Ángulo \phi','FontSize',TamanioFuente);
plot(t,fi,'*');
subplot(3,2,4);hold on;grid on; title('Velocidad ángulo \phi_p','FontSize',TamanioFuente);
plot(t,fi_p,'*');
subplot(3,2,1);hold on; grid on;title('Posición carro \delta','FontSize',TamanioFuente);hold on;
plot(t,d,'*');
subplot(3,2,3);hold on; grid on;title('Velocidad carro \delta_p','FontSize',TamanioFuente);hold on;
plot(t,d_p,'*');
subplot(3,1,3); grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
plot(t,u,'*');