%classic music
clear all
clc
tic
%参数设定
M=10;
doa=[-10 45 60]/180*pi;
P=length(doa);
f=1000;c=1500;
lambda=f/c;d=lambda/2;
snr=5;
N=400;
%阵列流型A
for i=1:P
    A(:,i)=exp(-j*2*pi*d*[0:M-1]'/lambda*sin(doa(i)));
end
%信源模型建立
for k=1:P
    S(k,:)=sqrt(10.^(snr/10))*(sind(1:N)+rand(1,N));
end
%接收信号模型建立
X=A*S;
%协方差矩阵特征值分解得到噪声子空间
R=X*X'/N;
[V,D]=eig(R);
[Y,I]=sort(diag(D));
Un=V(:,I(1:M-P));
%谱峰搜索部分
theta=-90:0.1:90;%线阵的搜索范围为-90~90度
 for i=1:length(theta)
   a_theta=exp(-j*(0:M-1)'*2*pi*d*sin(pi*theta(i)/180)/lambda);
   Pmusic(i)=1./abs((a_theta)'*Un*Un'*a_theta);
 end
plot(theta,10*log(Pmusic/max(Pmusic)),'r');
axis([-90 90 -90 10]);
xlabel('theta/degree');
ylabel('归一化空间谱/dB');
legend('Music Spectrum');
title('经典MUSIC估计');
grid on;
toc