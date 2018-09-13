clear all;
close all;
clc;
%����ʵ������
M = 4; N=5;
S=[(0:N-1)*M,(1:M-1)*N].';
LEN_S = length(S);

%����Դ�ź�
theta_bar = linspace(-60,60,25);  %ʵ��DOA�Ƕ�
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

%ʵ�ʽ����ź�
x_S = A_S * Source + noise_std * Noise;

%�㷨����
R_S1 = x_S * x_S' / SNAPSHOTS;    %���
R_S2 = x_S * x_S.' / SNAPSHOTS;    %�����
R_S3 = conj(x_S) * x_S' / SNAPSHOTS;    %�����
vecULA = getMaxULA(S);

%�ؽ�Z����
Z_vec = (rebuildZ([R_S1,R_S2,R_S3],vecULA,S)).';

%�ռ�ƽ������
DOF = 29;  %�������
m = length(Z_vec) - DOF + 1;  %ÿ��������Ԫ��
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
xlabel('�����/��');
ylabel('�׷�dB');
title('ǰ��ռ�ƽ��MUSIC����');

ylim=get(gca,'Ylim'); % ��ȡ��ǰͼ�ε�����ķ�Χ
hold on
plot([theta_bar;theta_bar],ylim,'b--'); % ����DOA�Ĳο�ֱ��
legend('FORWARD-SMOOTHNESS-MUSIC','REAL-DOA');
%ylim ���ڻ���y���ȡֵ��Χ
grid on;