clc
clear all
cpu=cputime;
new = zeros;
%Сигнал 'blocks'
ind = linspace(0,1,2^10); 
%Соотношение сигнал-шум 100
[x,noisyx] = wnoise(1,10,10);%Генерирование сигнала с аддитивным белым Гауссовым шумом
% noisyx = noisyx./max(noisyx);
figure(1)
subplot(2,1,1);
plot(ind,x)
grid on
v1=[0.0625 0.0625];
t1=[-10 20];
v2=[0.0625*2 0.0625*2];
t2=[-10 20];
v3=[0.0625*3 0.0625*3];
t3=[-10 20];
v4=[0.0625*4 0.0625*4];
t4=[-10 20];
v5=[0.0625*5 0.0625*5];
t5=[-10 20];
v6=[0.0625*6 0.0625*6];
t6=[-10 20];
v7=[0.0625*7 0.0625*7];
t7=[-10 20];
v8=[0.0625*8 0.0625*8];
t8=[-10 20];
v9=[0.0625*9 0.0625*9];
t9=[-10 20];
v10=[0.0625*10 0.0625*10];
t10=[-10 20];
v11=[0.0625*11 0.0625*11];
t11=[-10 20];
v12=[0.0625*12 0.0625*12];
t12=[-10 20];
v13=[0.0625*13 0.0625*13];
t13=[-10 20];
v14=[0.0625*14 0.0625*14];
t14=[-10 20];
v15=[0.0625*15 0.0625*15];
t15=[-10 20];
v16=[0.0625*16 0.0625*16];
t16=[-10 20];
subplot(2,1,2),plot(ind,noisyx,v1,t1,'black',v2,t2,'black',v3,t3,'black',v4,t4,'black',v5,t5,'black',v6,t6,'black',v7,t7,'black',v8,t8,'black',v9,t9,'black',v10,t10,'black',v11,t11,'black',v12,t12,'black',v13,t13,'black',v14,t14,'black',v15,t15,'black',v16,t16,'black')
grid on
hold on

%Поиск коэффициентов ДКП зашумленного сигнала
N=length(noisyx);
j=1;
w=1;
s=0;
DKP=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end

%Выполнение первого задания
DKPabs=abs(DKP);%Модуль ДКП коэффициентов
DKPabs(1)=[];
M = median(DKPabs);%Медиана
rezult1=M*1.483;%Результат первого задания
fprintf('Результат первого измерения = %f\n',rezult1);

%Выполнение второго задания
noisyx1=noisyx(1:512);
N=length(noisyx1);
j=1;
w=1;
s=0;
DKP1=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx1(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP1(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx2=noisyx(513:1024);
N=length(noisyx2);
j=1;
w=1;
s=0;
DKP2=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx2(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP2(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
DKPabs1=abs(DKP1);
DKPabs2=abs(DKP2);
DKPabs1(1)=[];
M1 = median(DKPabs1);
DKPabs2(1)=[];
M2 = median(DKPabs2);
M=(M1+M2)/2;
rezult1=M*1.483;%Результат
fprintf('Результат второго измерения = %f\n',rezult1);

%Выполнение третьего задания
noisyx1=noisyx(1:256);
N=length(noisyx1);
j=1;
w=1;
s=0;
DKP1=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx1(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP1(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx2=noisyx(257:512);
N=length(noisyx2);
j=1;
w=1;
s=0;
DKP2=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx2(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP2(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx3=noisyx(513:768);
N=length(noisyx3);
j=1;
w=1;
s=0;
DKP3=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx3(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP3(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx4=noisyx(769:end);
N=length(noisyx4);
j=1;
w=1;
s=0;
DKP4=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx4(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP4(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
DKPabs1=abs(DKP1);
DKPabs2=abs(DKP2);
DKPabs3=abs(DKP3);
DKPabs4=abs(DKP4);
Med=zeros;
DKPabs1(1)=[];
Med(1) = median(DKPabs1);
DKPabs2(1)=[];
Med(2) = median(DKPabs2);
DKPabs3(1)=[];
Med(3) = median(DKPabs3);
DKPabs4(1)=[];
Med(4) = median(DKPabs4);
M=mean(Med);%Среднее медиан
rezult1=M*1.483;%Результат
fprintf('Результат третьего измерения = %f\n',rezult1);

%Выполнение четвертого задания
noisyx1=noisyx(1:256);
N=length(noisyx1);
j=1;
w=1;
s=0;
DKP1=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx1(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP1(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx2=noisyx(257:512);
N=length(noisyx2);
j=1;
w=1;
s=0;
DKP2=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx2(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP2(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx3=noisyx(513:768);
N=length(noisyx3);
j=1;
w=1;
s=0;
DKP3=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx3(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP3(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx4=noisyx(769:end);
N=length(noisyx4);
j=1;
w=1;
s=0;
DKP4=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx4(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP4(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
DKPabs1=abs(DKP1);
DKPabs2=abs(DKP2);
DKPabs3=abs(DKP3);
DKPabs4=abs(DKP4);
Med=zeros;
DKPabs1(1)=[];
Med(1) = median(DKPabs1);
DKPabs2(1)=[];
Med(2) = median(DKPabs2);
DKPabs3(1)=[];
Med(3) = median(DKPabs3);
DKPabs4(1)=[];
Med(4) = median(DKPabs4);
M=median(Med);%Медиана медиан
rezult1=M*1.483;%Результат
fprintf('Результат четвертого измерения = %f\n',rezult1);

%Выполнение пятого задания
noisyx1=noisyx(1:128);
N=length(noisyx1);
j=1;
w=1;
s=0;
DKP1=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx1(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP1(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx2=noisyx(129:256);
N=length(noisyx2);
j=1;
w=1;
s=0;
DKP2=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx2(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP2(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx3=noisyx(257:384);
N=length(noisyx3);
j=1;
w=1;
s=0;
DKP3=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx3(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP3(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx4=noisyx(354:512);
N=length(noisyx4);
j=1;
w=1;
s=0;
DKP4=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx4(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP4(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end


noisyx5=noisyx(513:640);
N=length(noisyx5);
j=1;
w=1;
s=0;
DKP5=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx5(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP5(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx6=noisyx(641:768);
N=length(noisyx6);
j=1;
w=1;
s=0;
DKP6=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx6(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP6(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx7=noisyx(769:896);
N=length(noisyx7);
j=1;
w=1;
s=0;
DKP7=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx7(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP7(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx8=noisyx(897:end);
N=length(noisyx8);
j=1;
w=1;
s=0;
DKP8=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx8(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP8(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
DKPabs1=abs(DKP1);
DKPabs2=abs(DKP2);
DKPabs3=abs(DKP3);
DKPabs4=abs(DKP4);

DKPabs5=abs(DKP5);
DKPabs6=abs(DKP6);
DKPabs7=abs(DKP7);
DKPabs8=abs(DKP8);
Med=zeros;
DKPabs1(1)=[];
Med(1) = median(DKPabs1);
DKPabs2(1)=[];
Med(2) = median(DKPabs2);
DKPabs3(1)=[];
Med(3) = median(DKPabs3);
DKPabs4(1)=[];
Med(4) = median(DKPabs4);
DKPabs5(1)=[];
Med(5) = median(DKPabs5);
DKPabs6(1)=[];
Med(6) = median(DKPabs6);
DKPabs7(1)=[];
Med(7) = median(DKPabs7);
DKPabs8(1)=[];
Med(8) = median(DKPabs8);
M=mean(Med);
rezult1=M*1.483;%Результат
fprintf('Результат пятого измерения = %f\n',rezult1);

%Выполнение шестого задания
noisyx1=noisyx(1:128);
N=length(noisyx1);
j=1;
w=1;
s=0;
DKP1=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx1(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP1(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx2=noisyx(129:256);
N=length(noisyx2);
j=1;
w=1;
s=0;
DKP2=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx2(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP2(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx3=noisyx(257:384);
N=length(noisyx3);
j=1;
w=1;
s=0;
DKP3=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx3(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP3(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx4=noisyx(354:512);
N=length(noisyx4);
j=1;
w=1;
s=0;
DKP4=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx4(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP4(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end


noisyx5=noisyx(513:640);
N=length(noisyx5);
j=1;
w=1;
s=0;
DKP5=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx5(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP5(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx6=noisyx(641:768);
N=length(noisyx6);
j=1;
w=1;
s=0;
DKP6=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx6(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP6(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx7=noisyx(769:896);
N=length(noisyx7);
j=1;
w=1;
s=0;
DKP7=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx7(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP7(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx8=noisyx(897:end);
N=length(noisyx8);
j=1;
w=1;
s=0;
DKP8=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx8(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP8(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
DKPabs1=abs(DKP1);
DKPabs2=abs(DKP2);
DKPabs3=abs(DKP3);
DKPabs4=abs(DKP4);

DKPabs5=abs(DKP5);
DKPabs6=abs(DKP6);
DKPabs7=abs(DKP7);
DKPabs8=abs(DKP8);
Med=zeros;
DKPabs1(1)=[];
Med(1) = median(DKPabs1);
DKPabs2(1)=[];
Med(2) = median(DKPabs2);
DKPabs3(1)=[];
Med(3) = median(DKPabs3);
DKPabs4(1)=[];
Med(4) = median(DKPabs4);
DKPabs5(1)=[];
Med(5) = median(DKPabs5);
DKPabs6(1)=[];
Med(6) = median(DKPabs6);
DKPabs7(1)=[];
Med(7) = median(DKPabs7);
DKPabs8(1)=[];
Med(8) = median(DKPabs8);
M=median(Med);
rezult1=M*1.483;%Результат
fprintf('Результат шестого измерения = %f\n',rezult1);

%Выполнение седьмого задания
noisyx1=noisyx(1:64);
N=length(noisyx1);
j=1;
w=1;
s=0;
DKP1=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx1(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP1(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx2=noisyx(65:128);
N=length(noisyx2);
j=1;
w=1;
s=0;
DKP2=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx2(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP2(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx3=noisyx(129:192);
N=length(noisyx3);
j=1;
w=1;
s=0;
DKP3=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx3(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP3(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx4=noisyx(193:256);
N=length(noisyx4);
j=1;
w=1;
s=0;
DKP4=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx4(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP4(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end


noisyx5=noisyx(257:320);
N=length(noisyx5);
j=1;
w=1;
s=0;
DKP5=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx5(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP5(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx6=noisyx(321:384);
N=length(noisyx6);
j=1;
w=1;
s=0;
DKP6=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx6(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP6(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx7=noisyx(385:448);
N=length(noisyx7);
j=1;
w=1;
s=0;
DKP7=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx7(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP7(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx8=noisyx(449:512);
N=length(noisyx8);
j=1;
w=1;
s=0;
DKP8=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx8(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP8(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end




noisyx9=noisyx(513:576);
N=length(noisyx9);
j=1;
w=1;
s=0;
DKP9=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx9(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP9(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx10=noisyx(577:640);
N=length(noisyx10);
j=1;
w=1;
s=0;
DKP10=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx10(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP10(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx11=noisyx(641:704);
N=length(noisyx11);
j=1;
w=1;
s=0;
DKP11=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx11(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP11(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx12=noisyx(705:768);
N=length(noisyx12);
j=1;
w=1;
s=0;
DKP12=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx12(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP12(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end


noisyx13=noisyx(769:832);
N=length(noisyx13);
j=1;
w=1;
s=0;
DKP13=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx13(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP13(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx14=noisyx(833:896);
N=length(noisyx14);
j=1;
w=1;
s=0;
DKP14=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx14(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP14(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx15=noisyx(897:960);
N=length(noisyx15);
j=1;
w=1;
s=0;
DKP15=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx15(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP15(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx16=noisyx(961:end);
N=length(noisyx16);
j=1;
w=1;
s=0;
DKP16=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx16(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP16(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
DKPabs1=abs(DKP1);
DKPabs2=abs(DKP2);
DKPabs3=abs(DKP3);
DKPabs4=abs(DKP4);

DKPabs5=abs(DKP5);
DKPabs6=abs(DKP6);
DKPabs7=abs(DKP7);
DKPabs8=abs(DKP8);

DKPabs9=abs(DKP9);
DKPabs10=abs(DKP10);
DKPabs11=abs(DKP11);
DKPabs12=abs(DKP12);

DKPabs13=abs(DKP13);
DKPabs14=abs(DKP14);
DKPabs15=abs(DKP15);
DKPabs16=abs(DKP16);

Med=zeros;
DKPabs1(1)=[];
Med(1) = median(DKPabs1);
DKPabs2(1)=[];
Med(2) = median(DKPabs2);
DKPabs3(1)=[];
Med(3) = median(DKPabs3);
DKPabs4(1)=[];
Med(4) = median(DKPabs4);
DKPabs5(1)=[];
Med(5) = median(DKPabs5);
DKPabs6(1)=[];
Med(6) = median(DKPabs6);
DKPabs7(1)=[];
Med(7) = median(DKPabs7);
DKPabs8(1)=[];
Med(8) = median(DKPabs8);
DKPabs9(1)=[];
Med(9) = median(DKPabs9);
DKPabs10(1)=[];
Med(10) = median(DKPabs10);
DKPabs11(1)=[];
Med(11) = median(DKPabs11);
DKPabs12(1)=[];
Med(12) = median(DKPabs12);
DKPabs13(1)=[];
Med(13) = median(DKPabs13);
DKPabs14(1)=[];
Med(14) = median(DKPabs14);
DKPabs15(1)=[];
Med(15) = median(DKPabs15);
DKPabs16(1)=[];
Med(16) = median(DKPabs16);
M=mean(Med);
rezult1=M*1.483;%Результат
fprintf('Результат седьмого измерения = %f\n',rezult1);

%Выполнение восьмого задания
noisyx1=noisyx(1:64);
N=length(noisyx1);
j=1;
w=1;
s=0;
DKP1=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx1(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP1(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx2=noisyx(65:128);
N=length(noisyx2);
j=1;
w=1;
s=0;
DKP2=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx2(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP2(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx3=noisyx(129:192);
N=length(noisyx3);
j=1;
w=1;
s=0;
DKP3=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx3(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP3(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx4=noisyx(193:256);
N=length(noisyx4);
j=1;
w=1;
s=0;
DKP4=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx4(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP4(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end


noisyx5=noisyx(257:320);
N=length(noisyx5);
j=1;
w=1;
s=0;
DKP5=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx5(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP5(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx6=noisyx(321:384);
N=length(noisyx6);
j=1;
w=1;
s=0;
DKP6=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx6(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP6(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx7=noisyx(385:448);
N=length(noisyx7);
j=1;
w=1;
s=0;
DKP7=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx7(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP7(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx8=noisyx(449:512);
N=length(noisyx8);
j=1;
w=1;
s=0;
DKP8=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx8(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP8(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end




noisyx9=noisyx(513:576);
N=length(noisyx9);
j=1;
w=1;
s=0;
DKP9=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx9(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP9(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx10=noisyx(577:640);
N=length(noisyx10);
j=1;
w=1;
s=0;
DKP10=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx10(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP10(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx11=noisyx(641:704);
N=length(noisyx11);
j=1;
w=1;
s=0;
DKP11=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx11(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP11(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx12=noisyx(705:768);
N=length(noisyx12);
j=1;
w=1;
s=0;
DKP12=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx12(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP12(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end


noisyx13=noisyx(769:832);
N=length(noisyx13);
j=1;
w=1;
s=0;
DKP13=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx13(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP13(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx14=noisyx(833:896);
N=length(noisyx14);
j=1;
w=1;
s=0;
DKP14=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx14(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP14(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx15=noisyx(897:960);
N=length(noisyx15);
j=1;
w=1;
s=0;
DKP15=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx15(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP15(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
noisyx16=noisyx(961:end);
N=length(noisyx16);
j=1;
w=1;
s=0;
DKP16=zeros;
for k=0:1:N-1
    x1=sqrt((2/N));
    if k==0
        x1=1/(sqrt(N));
    end
    for n=0:1:N-1
        s=s+noisyx16(j)*cos(((2*n+1)*k*pi)/(2*N));
        j=j+1;
    end
    DKP16(w)=s*x1;
    w=w+1;
    j=1;
    s=0;
end
DKPabs1=abs(DKP1);
DKPabs2=abs(DKP2);
DKPabs3=abs(DKP3);
DKPabs4=abs(DKP4);

DKPabs5=abs(DKP5);
DKPabs6=abs(DKP6);
DKPabs7=abs(DKP7);
DKPabs8=abs(DKP8);

DKPabs9=abs(DKP9);
DKPabs10=abs(DKP10);
DKPabs11=abs(DKP11);
DKPabs12=abs(DKP12);

DKPabs13=abs(DKP13);
DKPabs14=abs(DKP14);
DKPabs15=abs(DKP15);
DKPabs16=abs(DKP16);

Med=zeros;
DKPabs1(1)=[];

figure(2);
% DKP2(1)=[];
% col=hist(DKP16);
% max=max(col);
% max1=-20:20;
% histogram(DKP2,max1)

Med(1) = median(DKPabs1);
DKPabs2(1)=[];
Med(2) = median(DKPabs2);
DKPabs3(1)=[];
Med(3) = median(DKPabs3);
DKPabs4(1)=[];
Med(4) = median(DKPabs4);
DKPabs5(1)=[];
Med(5) = median(DKPabs5);
DKPabs6(1)=[];
Med(6) = median(DKPabs6);
DKPabs7(1)=[];
Med(7) = median(DKPabs7);
DKPabs8(1)=[];
Med(8) = median(DKPabs8);
DKPabs9(1)=[];
Med(9) = median(DKPabs9);
DKPabs10(1)=[];
Med(10) = median(DKPabs10);
DKPabs11(1)=[];
Med(11) = median(DKPabs11);
DKPabs12(1)=[];
Med(12) = median(DKPabs12);
DKPabs13(1)=[];
Med(13) = median(DKPabs13);
DKPabs14(1)=[];
Med(14) = median(DKPabs14);
DKPabs15(1)=[];
Med(15) = median(DKPabs15);
DKPabs16(1)=[];
Med(16) = median(DKPabs16);
M=median(Med);
rezult1=M*1.483;%Результат
Medblock=1.483*Med;
fprintf('Результат восьмого измерения = %f\n',rezult1);

bar(Medblock)
grid on
cpu1=cputime;
disp(cpu1-cpu)
