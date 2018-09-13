clear all;
close all;
clc;
%设置实际阵列
M = 4; N=5;
S=[(0:N-1)*M,(1:M-1)*N].';
LEN_S = length(S);

%生成源信号
theta_bar = linspace(-60,60,25);  %实际DOA角度
DD = length(theta_bar);
SNRdB = 0;
SNAPSHOTS = 2000;
d_lumbda = 1/2;
A_S = exp(2i * pi * d_lumbda * S * sind(theta_bar));
for jj = 1:DD
Source(jj,:) = sin(2 * pi * randn(1, SNAPSHOTS));
end
Noise = (randn(LEN_S, SNAPSHOTS) + 1i * randn(LEN_S, SNAPSHOTS)) / sqrt(2);
noise_std = 10^(-SNRdB/20);

%实际接收信号
x_S = A_S * Source + noise_std * Noise;

%算法处理
R_S1 = x_S * x_S' / SNAPSHOTS;    %相减
R_S2 = x_S * x_S.' / SNAPSHOTS;    %正相加
R_S3 = conj(x_S) * x_S' / SNAPSHOTS;    %负相加
vecULA = getMaxULA(S);

%重建Z向量
Z_vec = (rebuildZ([R_S1,R_S2,R_S3],vecULA,S)).';

%空间平滑处理
DOF = 29;  %子阵个数
m = length(Z_vec) - DOF + 1;  %每个子阵阵元数
R_zz = zeros(m,m);
for j = 1:DOF
    Z_i  = Z_vec(j:j+m-1);
    R_zz = R_zz + Z_i*Z_i';
end
R_zz = R_zz/DOF;

%MUSIC algorithm
[U,Ss,V]=svd(R_zz);
Un=U(:,DD+1:m);

theta=-60:0.1:60;
for ii=1:length(theta)
    a=exp(2i * pi *  d_lumbda * (0:m-1)' * sind(theta(ii)));
    Psmoothmusic(ii)=1./abs(a'*Un*Un'*a);
end

plot(theta,10*log(Psmoothmusic)','-r');
xlabel('入射角/度');
ylabel('谱峰dB');
title('前向空间平滑MUSIC估计');

ylim=get(gca,'Ylim'); % 获取当前图形的纵轴的范围
hold on
plot([theta_bar;theta_bar],ylim,'b--'); % 绘制DOA的参考直线
legend('FORWARD-SMOOTHNESS-MUSIC','REAL-DOA');
%ylim 用于绘制y轴的取值范围
grid on;