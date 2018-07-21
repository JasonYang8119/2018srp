% Sample codes for 
%   1) the coarray MUSIC algorithm [1] and 
%   2) coarray interpolation [2].
%
%	Written by Chun-Lin Liu
%	E-mail : cl.liu@caltech.edu
%	Project website:
%	http://systems.caltech.edu/dsp/students/clliu/coarray.html
%	Last revised on July 13, 2017
%
%   If you use this script, please cite reference [1] or [2].
%
%	Reference
%   [1] C.-L. Liu and P. P. Vaidyanathan, â€œRemarks on the Spatial Smoothing Step in Coarray MUSIC,â€?IEEE Signal Processing Letters, vol. 22, no. 9, pp. 1438-1442, Sep. 2015. 
%   [2] C.-L. Liu, P. P. Vaidyanathan and P. Pal, â€œCoprime Coarray Interpolation for DOA Estimation via Nuclear Norm Minimization,â€?in Proc. of 2016 IEEE International Symposium on Circuits and Systems (ISCAS 2016), pp. 2639-2642, Montreal, Canada, May 2016. 
% 


clear; clc; close all;


%% ========== Array configurations ==========
% ========== Coprime array ==========
M = 4; N = 5;
% Sensor locations
S=[(0:N-1)*M,(1:M-1)*N].';
LEN_S = length(S);
%% ========== Generate array outputs x_S ==========
% True normalized DOA %sin(-60¡ã~60¡ã)
theta_bar = linspace(-60, 60, 25);
% Number of sources
DD = length(theta_bar);
% Number of snapshots
SNAPSHOTS = 2000;
% SNR in dB
SNRdB = 0;

% Array manifold matrix V_S
A_S = exp(2i * pi * S * sin(theta_bar*pi/180));

Source = sin ( pi* randn(DD, SNAPSHOTS));
Noise = randn(LEN_S, SNAPSHOTS);

noise_std = 10^(-SNRdB/20);

% Array output
x_S = A_S * Source + noise_std * Noise;


%% ========== Sample covariance matrix R_S ==========
R_S1 = x_S * x_S' / SNAPSHOTS;
R_S2 = x_S * x_S.' / SNAPSHOTS;
R_S3 = conj(x_S) * (conj(x_S))' / SNAPSHOTS;

%% ========== Convert R_S to sample autocorrelation vectors on the difference&sum coarray ==========
%virtual co-array
[D,~] = weight_function(S,'SD');
%Z=vec(R_S)
R_S_all = [R_S2 R_S1 R_S3];
S = unique(S);
[n1, n2] = ndgrid(S);
n1_n2_sd = [n1+n2 n1-n2 -n1-n2];
LEN_SD = length(D);
Z_SD = zeros( LEN_SD, 1 );
for mm = 1 : LEN_SD
    % Select samples with the same lag
    data = R_S_all( n1_n2_sd == D(mm) );
    % ====================
    % Average over all the samples with the same lag.
    Z_SD(mm) = mean( data );
end
%% ======== spacial smoothing technique ========
%number of sub-array and length of each sub_array
N_sub = (LEN_SD+1)/2;
sub_LEN = LEN_SD + 1 - N_sub;
%temple sub_array
Zi_SD = zeros(sub_LEN,1);
%DOF equals number of sub_array
R_Zi = zeros(sub_LEN,sub_LEN,N_sub);
for ii = 1:N_sub 
    %divide into multiple overlapping sub-arrays
    Zi_SD= Z_SD(ii:(ii+sub_LEN-1));
    R_Zi(:,:,ii) = Zi_SD * Zi_SD.';
end
%get mean of each sub_array
R_cov = mean(R_Zi,3);
%% ========== Coarray MUSIC algorithm ==========
[EV, EW] = eig(R_cov);
% Sort by the absolute values of eigenvalues
ew = diag(EW);%ew is eigenvalues column vector EV is eigenvalues matrix
[~, II] = sort(abs(ew), 'descend');
% Sort eigenvectors
% U_all = [Us, Un]
% Us and Un are the signal and noise subspace, respectively.
U_all = EV(:, II);
% Check whether DD exceeds the dimension
if (DD > sub_LEN && DD > N_sub)
    disp('Warning, the number of sources exceeds the limit of sub-coarry. Choose only one vector as noise subspace.');
    Un = U_all(:, end);
else
    Un = U_all(:, DD+1:end); % Noise subspace
end
%DOA spectrum search
iii = 1;
for search_i=-60:0.01:60
a_vS = exp(2i * pi *(-28:0).'* sin(search_i*pi/180));
P_msc(iii) = 1/((a_vS'*Un)*(Un'*a_vS));
iii = iii + 1;
end
P_msc = P_msc/max(P_msc);
%% ========== MUSIC Spectra ==========
plot(-60:0.01:60,P_msc)
