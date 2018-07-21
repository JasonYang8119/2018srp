% Sensor locations
S = [-4:4].';
% Number of sensors
LEN_S = length(S);
%test MUSIC
%True normalized DOA
bar = linspace(-90, 90, 5).';
theta_bar = sin(pi*bar/180);
% Number of sources
DD = length(theta_bar);
% Number of snapshots
SNAPSHOTS = 500;
% SNR in dB
SNRdB = 0;
% Array manifold matrix V_S
V_S = exp(1i * pi * S * theta_bar.');
Source = (randn(DD, SNAPSHOTS) + 1i * randn(DD, SNAPSHOTS)) / sqrt(2);
Noise = (randn(LEN_S, SNAPSHOTS) + 1i * randn(LEN_S, SNAPSHOTS)) / sqrt(2);
noise_std = 10^(-SNRdB/20);
% Array output
x_S = V_S * Source + noise_std * Noise;
R_S = x_S * x_S' / SNAPSHOTS;
[P,theta] = clscl_MUSIC(R_S,S,DD,bar);
plot(theta,P)

% [D,w] = weight_function(S,'SD');
% [n1, n2] = ndgrid(S);
% n1_n2_mat = n1 - n2;
% n1_n2_vec = n1_n2_mat(:);
% n1_n2_pls = n1 + n2;
% n1_n2_plsv = n1_n2_pls(:);
% N1_n = 2*[(0:N-1)*M].';
% N2_n = 2*[(0:M-1)*N].';
% n1_n2_vecp = [-n1_n2_plsv.' -N1_n.' -N2_n.' n1_n2_vec.' n1_n2_plsv.' N1_n.' N2_n.'].';
% n1_n2_pl = unique(n1_n2_vecp, 'sorted');
% max_pl = max(n1_n2_pl);
% max_vec = -max_pl:max_pl;
% fin_hole = setdiff(max_vec,n1_n2_pl);
% min_hole = min(abs(fin_hole));
% [U, ~] = weight_function( S, 'U' );
% x_U = zeros( length(U), 1 );
% for mm=1:length(U)
%  data = R_S(n1 - n2 == U(mm));
%  x_U(mm) = mean( data );
%  len = (n1 - n2 == U(mm));
% end

% Number of snapshots
% SNAPSHOTS = 2000;
% 
% Source = sin(pi * randn(25, SNAPSHOTS));
% Input_S = Source';
% data = mean(Input_S);

% test how to get a matrix
% A = [1,2,3;4,5,6;7,8,9;10,11,12];
% N1=ones(4,3);
% N2=[0,1,1;1,0,1;1,1,0;0,0,0];
% B=A(N1-N2==1)

% %SD = negtive_sum + positive_sum +difference
% %difference have beed done: 
% n1_n2_mat = n1 - n2;
% n1_n2_vec = n1_n2_mat(:);
% %positive_sum
% n1_n2_psum = n1 + n2;
% %negtive_sum
% n1_n2_nsum = -n1 - n2;
% n1_n2_mat = [n1_n2_psum n1_n2_mat n1_n2_nsum];
% %new vector
% n1_n2_vec = n1_n2_mat(:);
% D_MAX = max(abs(n1_n2_vec));
% D_full = (-D_MAX : D_MAX).';
% D_set = unique(n1_n2_vec, 'sorted');
% w_set = histc(n1_n2_vec, D_set);
% H = setdiff(D_full, D_set); % Holes
% if (isempty(H))
%   % D is hole-free
%   D = D_set; w = w_set;
% else
% N_max = min(abs(H)) - 1;
% D = (-N_max : N_max).';
% w = w_set(D_set >= -N_max & D_set <= N_max);
% end