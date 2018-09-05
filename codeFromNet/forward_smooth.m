%FORWARD_SMOOTHNESS_MUSIC ALOGRITHM VS CLASSIC MUSIC
%DOA ESTIMATION BY FORWARD_SMOOTHNESS_MUSIC
clear all;
close all;
clc;

M=8;%阵元数
N=1024; %信号长度
snapshot_number=N;%快拍数
doa=[45,-60];%信号波达方向
P=length(doa);%信元数
f0=1000;%信号频率
fs=200;Ts=1/fs;
c=1500;%水中声速1500m/s
l=c/f0;%信号波长  
d=0.5*l;%阵元间距
m=6;%每个子阵阵元数
L=M-m+1;%相互交错的子阵数
snr=0;%信噪比

A=[exp(-j*(0:M-1)*d*2*pi*sin(doa(1)*pi/180)/l);exp(-j*(0:M-1)*d*2*pi*sin(doa(2)*pi/180)/l)].';%阵列流型

t=(0:N-1)*Ts;
s1=sqrt(10.^(snr/10))*exp(j*(2*pi*f0*t+pi/1));
s2=sqrt(10.^(snr/10))*exp(j*(2*pi*f0*t+pi/2));
s=[s1;s2];%仿真信号,显然这是一个相干信源
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(M,N)+j*randn(M,N));%加了高斯白噪声后的阵列接收信号
R=x*x'/N;
Rf=zeros(m,m);
for i=1:L
    y=x(i:i+m-1,:);
    Rf=Rf+y*y'/N;
end
Rf=Rf/L;

[U,S,V]=svd(Rf);
Un=U(:,P+1:m);

[V,D]=eig(R);
UN=V(:,1:M-P);


theta=-90:0.1:90;%线阵的搜索范围为-90~90度
 for ii=1:length(theta)
   a=exp(-j*(0:m-1)'*2*pi*d*sin(pi*theta(ii)/180)/l);
   Psmoothmusic(ii)=1./abs(a'*Un*Un'*a);
 end
 
 for kk=1:length(theta)
   a=exp(-j*(0:M-1)'*2*pi*d*sin(pi*theta(kk)/180)/l);
   Pmusic(kk)=1./abs(a'*UN*UN'*a);
 end

plot(theta,10*log(Psmoothmusic),'R-',theta,10*log(Pmusic),'B-');
%axis([-90 90 -90 90]);
xlabel('入射角/度');
ylabel('谱峰dB');
legend('FORWARD-SMOOTHNESS-MUSIC','CLASSIC MUSIC');
title('前向空间平滑MUSIC估计与经典MUSIC的对比');
grid on;