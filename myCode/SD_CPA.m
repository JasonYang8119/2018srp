clear all;
close all;
clc;
%Set sensor
M = 4; N=5;
S=[(0:N-1)*M,(1:M-1)*N].';  %column-vector,true field array set
LEN_S = length(S);

%generate sinusoidal sources
theta_bar = linspace(-60,60,5);  %row-vector,true DOA angel
% Number of sources
DD = length(theta_bar);
SNRdB = 0;
SNAPSHOTS = 2000;
A_S = exp(2i * pi * S * sind(theta_bar));

Source = randn(DD, SNAPSHOTS);
Noise = (randn(LEN_S, SNAPSHOTS) + 1i * randn(LEN_S, SNAPSHOTS)) / sqrt(2);
noise_std = 10^(-SNRdB/20);

% Array output
x_S = A_S * Source + noise_std * Noise;

%autocorrelation 
R_S1 = x_S * x_S' / SNAPSHOTS;    %minus
R_S2 = x_S * x_S.' / SNAPSHOTS;    %add_positive
R_S3 = conj(x_S) * x_S' / SNAPSHOTS;    %add_negtive
vecULA = getMaxULA(S);

%Rebuild Z
Z_vec = (rebuildZ([R_S1,R_S2,R_S3],vecULA,S)).';

%spatial smoothing technique
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
   a=exp(2i * pi * (0:m-1)' * sind(theta(ii)));
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