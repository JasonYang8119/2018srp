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
%   [1] C.-L. Liu and P. P. Vaidyanathan, “Remarks on the Spatial Smoothing Step in Coarray MUSIC,” IEEE Signal Processing Letters, vol. 22, no. 9, pp. 1438-1442, Sep. 2015. 
%   [2] C.-L. Liu, P. P. Vaidyanathan and P. Pal, “Coprime Coarray Interpolation for DOA Estimation via Nuclear Norm Minimization,” in Proc. of 2016 IEEE International Symposium on Circuits and Systems (ISCAS 2016), pp. 2639-2642, Montreal, Canada, May 2016. 
% 


clear; clc; close all;


%% ========== Array configurations ==========
% % % ========== Nested array ==========
% % N1 = 3; N2 = 3;
% % % Sensor locations
% % S = [1:N1, (1:N2)*(N1+1)].';
% % % ==================================

% ========== Coprime array ==========
M = 5; N = 7;
% Sensor locations
S = [(0:N-1)*M, (1:2*M-1)*N].';
% ===================================

% Number of sensors
LEN_S = length(S);

%% ========== Generate array outputs x_S ==========
% True normalized DOA
theta_bar = linspace(-0.49, 0.49, 35).';
% Number of sources
DD = length(theta_bar);
% Number of snapshots
SNAPSHOTS = 500;
% SNR in dB
SNRdB = 0;

% Array manifold matrix V_S
V_S = exp(2i * pi * S * theta_bar.');

Source = (randn(DD, SNAPSHOTS) + 1i * randn(DD, SNAPSHOTS)) / sqrt(2);
Noise = (randn(LEN_S, SNAPSHOTS) + 1i * randn(LEN_S, SNAPSHOTS)) / sqrt(2);

noise_std = 10^(-SNRdB/20);

% Array output
x_S = V_S * Source + noise_std * Noise;


%% ========== Sample covariance matrix R_S ==========
R_S = x_S * x_S' / SNAPSHOTS;

%% ========== Convert R_S to sample autocorrelation vectors on the difference coarray ==========

% Compute the sample autocorrelation vector x_{S_{diff}^{ULA}} in [1]
[U, x_U] = sample_covariance_to_difference_coarray( S, R_S, 'x_U' );

% Compute the sample autocorrelation vector x_V after coarray interpolation
% in [2]
% Note: cvx package (http://cvxr.com/cvx/) is required here.
[V, x_V] = sample_covariance_to_difference_coarray( S, R_S, 'x_V' );


%% ========== Coarray MUSIC algorithm ==========

% Select central ULA segment [1]
[ theta_bar_grid_U, P_U, theta_bar_est_U ] = coarray_MUSIC( x_U, DD );

% Coarray interpolation [2]
[ theta_bar_grid_V, P_V, theta_bar_est_V ] = coarray_MUSIC( x_V, DD );


%% ========== MUSIC Spectra ==========

% Figure 2 in [1]
figure;
semilogy(theta_bar_grid_U, P_U)
xlabel('$\bar{\theta}$', 'interpreter', 'latex')
ylabel('Normalized MUSIC Spectrum', 'interpreter', 'latex')
set(gca, 'XTick', theta_bar);
grid on;
title(['RMSE = ', num2str( norm(theta_bar_est_U - theta_bar) / sqrt(DD) )])

% Similar to the results in [2]
figure;
semilogy(theta_bar_grid_V, P_V)
xlabel('$\bar{\theta}$', 'interpreter', 'latex')
ylabel('Normalized MUSIC Spectrum', 'interpreter', 'latex')
set(gca, 'XTick', theta_bar);
grid on;
title(['RMSE = ', num2str( norm(theta_bar_est_V - theta_bar) / sqrt(DD) )])

