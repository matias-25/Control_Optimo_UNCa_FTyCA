clear all;%close all; 
tic;
%% Constantes del sistema
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5; 

%Condiciones iniciales
d(1)=0;
d_p(1)=0;
% fi(1)=0.5;
fi(1)=0.1;
fi_p(1)=0;
u(1)=0; %accion de control 
x= [d(1);d_p(1);fi(1);fi_p(1)];
x0=x;

%Versión linealizada en el equilibrio inestable. Sontag Pp 104. 
Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0;0 0 0 1;0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_Bc=[0; 1/M; 0; -1/(long*M)];
Mat_Cc=[1 0 0 0; 0 0 1 0];

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
Q=diag([1e2 1e1 1e6 1e1]);R=1e1;%Matrices de diseño del controlador DLQG
% Q=diag([1e-6 1e-6 1e-3 1e-6]);R=1e6;
Kdlqr = dlqr(Mat_Ad,Mat_Bd,Q,R); %ganacia del controlador
disp('Polos de dlqr en:')
Ea=eig(Mat_Ad-Mat_Bd*Kdlqr)

%% Obsevador de Luenberger
%Cálculo del Observador--------------------------------------------------- 
A_o=Mat_Ad'; 
B_o=Mat_Cd'; 
C_o=Mat_Bd'; 
Qo=diag([1e0 1e0 1e0 1e0]);Ro=diag([1e-1, 1e-1]); 
% Qo=diag([1e0 1e0 1e0 1e0]);Ro=diag([1e0, 1e0]); 
Ko= (dlqr(A_o,B_o,Qo,Ro))'; %ganancia del observador
%%
forzador=1; %limit la accion de accion de control
Qcomp=eye(4);%Para comparar el desempeño de los controladores
%% Monte Carlo
% Consigna 0, 0.01, 0.02, 0.05 y 0.1
sQ=0.1; %Para F
sR=0.1; %Para G. Covarianza del ruido de medicion sigma=sqrt(sR)
F_=sQ*eye(4); %Covarianza del ruido de estado Sigma=sqrt(sQ)
G_=sR;
S=Q;
P=S; %condición inicial de P
kmax=2000;
Realizaciones=50; %Cantidad de realizaciones para el Monte Carlo.
Kx=zeros(kmax,4);
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
    Kx(hi,:)=inv(R+Ba'*P*Ba)*Ba'*P*Aa;
    Kv(hi,:)=inv(R+Ba'*P*Ba)*Ba'*P*F_;
    %Ea(:,hi)=eig(Aa-Ba*Kx(hi,:));
end

%% _____________ESTIMADOR KALMAN______________
P_Kalman=F_*F_';
P11=zeros(1,5000);P22=P11;P33=P11;P44=P11;
for h_k=1:5000
    P_Kalman_=Mat_Ad*P_Kalman*Mat_Ad'+(F_*F_');
    K_Kalman=P_Kalman_*Mat_Cd'/(Mat_Cd*P_Kalman_*Mat_Cd'+(G_*G_')); %Ganancia de Kalman
    P_Kalman=(eye(4)-K_Kalman*Mat_Cd)*P_Kalman_;
    P11(h_k)=P_Kalman(1,1);
    P22(h_k)=P_Kalman(2,2);
    P33(h_k)=P_Kalman(3,3);
    P44(h_k)=P_Kalman(4,4);
end
% figure(4);semilogy(P11);hold on;semilogy(P22,'g');semilogy(P33,'c');semilogy(P44,'r');title('Evolución de P_1_1,P_2_2,P_3_3 y P_4_4.')
% xlabel('Iteraciones');
disp('Polos de Kalman en:')
EK=abs(eig(Mat_Ad-K_Kalman*Mat_Cd))
%%
for trial=1:Realizaciones %Empieza el Monte Carlo
    v=randn(4,kmax);%Señales aleatorios de media nula y varianza unidad.
    w=randn(2,kmax);
    x=x0+F_*v(:,1);
    x_hat=[0;0;0;0];% variables estado observador 
    x_hat_=[0;0;0;0];
    d(trial,1)=x(1);
    d_p(trial,1)=x(2);
    fi(trial,1)=x(3);
    fi_p(trial,1)=x(4);
    
    for ki=1:kmax-1
        Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*Qcomp*x + u(trial,ki)'*R*u(trial,ki));
        Ys=Mat_Cd*x+G_*w(:,ki);
        y_sal(trial,ki)=Ys(1);
        x_hat=x_hat_+K_Kalman*(Ys -Mat_Cd*x_hat_);
%         u(trial,ki)=-Kx(1,:)*x_hat; %LQR
%         u(trial,ki)=-Kx(1,:) *x; %LQR
        uaux=-Kx(1,:) *x_hat; %LQR
        if abs(uaux)>forzador
            u(trial,ki)=forzador*sign(uaux);
        else
            u(trial,ki)=uaux;
        end 
        x=modelo_avion(Ts,x,u(trial,ki))+F_*v(:,ki+1);
        x_hat_=Mat_Ad*x_hat+Mat_Bd*u(trial,ki);
        
        d(trial,ki+1)=x(1);
        d_p(trial,ki+1)=x(2);
        fi(trial,ki+1)=x(3);
        fi_p(trial,ki+1)=x(4);
    end
    Jn_(trial,ki+1)=Jn_(trial,ki+1)+x'*S*x;
%     u(trial,ki+1)=-Kx(1,:)*[x_hat]-Kv(1,:)*[v(:,ki)];
%         u(trial,ki)=-Kx(1,:) *x; %LQR
        uaux=-Kx(1,:) *x_hat; %LQR
        if abs(uaux)>forzador
            u(trial,ki)=forzador*sign(uaux);
        else
            u(trial,ki)=uaux;
        end 
end

Jn=mean(Jn_);disp(['El valor de costo es Jn(end)=' num2str(Jn(end)) '.fi(1)=' num2str(fi(1)) '.']);
%% Graficos
t=t*Ts;
TamanioFuente=14;
figure;
subplot(3,2,2);hold on;grid on;title('Ángulo \phi','FontSize',TamanioFuente);
plot(t,mean(fi),'b');plot(t,mean(fi)+.5*sqrt(var(fi)),'r');plot(t,mean(fi)-.5*sqrt(var(fi)),'r');
subplot(3,2,4);hold on;grid on; title('Velocidad ángulo \phi_p','FontSize',TamanioFuente);
plot(t,mean(fi_p),'b');plot(t,mean(fi_p)+.5*sqrt(var(fi_p)),'r');plot(t,mean(fi_p)-.5*sqrt(var(fi_p)),'r');
subplot(3,2,1);hold on; grid on;title('Posición carro \delta','FontSize',TamanioFuente);hold on;
plot(t,mean(d),'b');plot(t,mean(d)+.5*sqrt(var(d)),'r');plot(t,mean(d)-.5*sqrt(var(d)),'r');
subplot(3,2,3);hold on; grid on;title('Velocidad carro \delta_p','FontSize',TamanioFuente);hold on;
plot(t,mean(d_p),'b');plot(t,mean(d_p)+.5*sqrt(var(d_p)),'r');plot(t,mean(d_p)-.5*sqrt(var(d_p)),'r');
subplot(3,1,3); grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
plot(t,mean(u),'b');plot(t,mean(u)+.5*sqrt(var(u)),'r');plot(t,mean(u)-.5*sqrt(var(u)),'r');
 disp(['tiempo de ejecucion :' num2str(toc) ' [seg].']);
% figure;hold on;
% subplot(2,2,1);hold on; plot(mean(fi),mean(fi_p),color); grid on;xlabel('Ángulo','FontSize',TamanioFuente);ylabel('Velocidad angular','FontSize',TamanioFuente);hold on;
% subplot(2,2,2);hold on; plot(d,d_p,color); grid on;xlabel('Posición carro','FontSize',TamanioFuente);ylabel('Velocidad carro','FontSize',TamanioFuente);hold on;
% subplot(2,2,3);hold on;
% plot(t,Jn,color);plot(t,Jmin*ones(size(t)),colorc);ylabel('Acumulación de costo','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);
% subplot(2,2,4);hold on; plot(abs(Ea)');ylabel('Polos de lazo cerrado','FontSize',TamanioFuente);xlabel('Etapas de iteración','FontSize',TamanioFuente);
