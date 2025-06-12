clear all;%close all; 
%% Constantes del sistema
% w= , a= , b= , c=velocidad del avion
w=2;a=0.05; b=5;c=80;

%Condiciones iniciales
alfa(1)=0;
fhi(1)=0;
fhi_p(1)=0;
h(1)=100;
u(1)=0; %accion de control 
x= [-alfa(1);-fhi(1);-fhi_p(1);-h(1)];
x0=-x;

%Versión linealizada del avion x=[alfa;fi;fi_p;h]
%alfa= direccion de vuelo, fi= angulo cabeceo (piloto) , 
%fi_p = velocidad angulo cabeceo, h= altura
Mat_Ac=[-a a 0 0;0 0 1 0; w^2 -w^2 0 0; c 0 0 0];
Mat_Bc=[0; 0; b*w^2; 0];
Mat_Cc=[0 0 0 1; %mido la altura
        0 1 0 0];   %mido angulo cabeceo

%% Sistema continuo
sys_c=ss(Mat_Ac,Mat_Bc,Mat_Cc,0); 

%% Sistema discreto
fs=500; %frecuencia de muestreo
Ts=1/fs;

sys_d=c2d(sys_c,Ts,'zoh'); 
Mat_Ad=sys_d.a;
Mat_Bd=sys_d.b;
Mat_Cd=sys_d.c;

%% controlador DLQG
Q=diag([1e2 1e2 1e0 5e-4]);R=1e0;%Matrices de diseño del controlador DLQG
Kx = dlqr(Mat_Ad,Mat_Bd,Q,R); %ganacia del controlador

%% Obsevador de Luenberger
%Cálculo del Observador--------------------------------------------------- 
A_o=Mat_Ad'; 
B_o=Mat_Cd'; 
C_o=Mat_Bd'; 
Qo=diag([1e-3 1e-3 1e-2 1e0]);Ro=diag([1e0, 1e0]); 
Ko= dlqr(A_o,B_o,Qo,Ro); %ganancia del observador
%%
Qcomp=eye(4);%Para comparar el desempeño de los controladores
%% Monte Carlo
% Consigna 0, 0.01, 0.02, 0.05 y 0.1
sQ=0.01; %Para F
sR=0.01; %Para G. Covarianza del ruido de medicion sigma=sqrt(sR)
F_=sQ*eye(4); %Covarianza del ruido de estado Sigma=sqrt(sQ)
G_=sR;
S=Q;
P=S; %condición inicial de P
kmax=20000;
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
randn('state',100);
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
    alfa(trial,1)=x(1);
    fhi(trial,1)=x(2);
    fhi_p(trial,1)=x(3);
    h(trial,1)=x(4);
    
    for ki=1:kmax-1
        %u(trial,ki)=-Kx(1,:)*x; %DLQR
        %Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*Q*x + u(trial,ki)'*R*u(trial,ki));
        u(trial,ki)=-Kx(1,:)*x_hat-Kv(1,:)*v(:,ki); %LQG Ec 6-45,con observador
        Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*Qcomp*x + u(trial,ki)'*R*u(trial,ki));%Ec 11-38 con observador
        %y_sal(trial,ki)=Mat_Cd*x+G_*w(ki);
        Ys=Mat_Cd*x + G_*w(:,ki);
        y_sal(trial,ki)=Ys(1);
        x=modelo_avion(Ts,x,u(trial,ki))+F_*v(:,ki+1);
        x_hat= Mat_Ad*x_hat + Mat_Bd*u(trial,ki) + Ko'*(Ys-Mat_Cd*x_hat);% Obsevador de Luenberger
        alfa(trial,ki+1)=x(1);
        fhi(trial,ki+1)=x(2);
        fhi_p(trial,ki+1)=x(3);
        h(trial,ki+1)=x(4);
    end
    Jn_(trial,ki+1)=Jn_(trial,ki+1)+x'*S*x;
    u(trial,ki+1)=-Kx(1,:)*[x]-Kv(1,:)*[v(:,ki)];
end

Jn=mean(Jn_);disp(['El valor de costo es Jn(end)=' num2str(Jn(end)) '.fi(1)=' num2str(fhi_p(1)) '.']);
%% Graficos
t=t*Ts;
TamanioFuente=14;
figure;
subplot(3,2,1);hold on;grid on; title('Altura avion h','FontSize',TamanioFuente);
plot(t,h,'*');
subplot(3,2,2);hold on; grid on;title('Ángulo cabeceo \phi','FontSize',TamanioFuente);hold on;
plot(t,fhi,'*');
subplot(3,2,4);hold on;grid on;title('Velocidad Ángulo cabeceo\phi_p','FontSize',TamanioFuente);
plot(t,fhi_p,'*');
subplot(3,2,3);hold on; grid on;title('Ángulo de vuelo \alpha','FontSize',TamanioFuente);hold on;
plot(t,alfa,'*');
subplot(3,1,3); grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
plot(t,u,'*');