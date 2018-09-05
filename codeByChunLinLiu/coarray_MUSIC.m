function varargout = coarray_MUSIC( x_U, DD )
%COARRAY_MUSIC: MUSIC spectrum using sample autocorrelation vector x_U
%==============================================================================
%
%	Written by Chun-Lin Liu
%	E-mail: cl.liu@caltech.edu
%	Project website:
%	http://systems.caltech.edu/dsp/students/clliu/coarray.html
%	Last revised on July 13, 2017
%
%   If you use this function, please cite reference [1].
%
%   Input:
%
%       1) x_U: Sample autocorrelation vector on the central ULA segment.
%           This is x_{S_{diff}^{ULA}} in [1].
%           This argument can be x_V after coarray interpolation [2].
%       2) DD: Number of sources
%
%   Output:
%
%       theta_bar_est = coarray_MUSIC( ... )
%       [theta_bar_grid, P] = coarray_MUSIC( ... )
%       [theta_bar_grid, P, theta_bar_est] = coarray_MUSIC( ... )
%
%       1) theta_bar_grid: A dense grid for the MUSIC spectrum.
%       2) P: normalized MUSIC spectrum with maximum 1.
%       3) theta_bar_est: The estimated normalized DOA using the root MUSIC
%       algorithm.
%
%   Reference
%   [1] C.-L. Liu and P. P. Vaidyanathan, “Remarks on the Spatial Smoothing Step in Coarray MUSIC,” IEEE Signal Processing Letters, vol. 22, no. 9, pp. 1438-1442, Sep. 2015.
%   [2] C.-L. Liu, P. P. Vaidyanathan and P. Pal, “Coprime Coarray Interpolation for DOA Estimation via Nuclear Norm Minimization,” in Proc. of 2016 IEEE International Symposium on Circuits and Systems (ISCAS 2016), pp. 2639-2642, Montreal, Canada, May 2016.
%
%==============================================================================


% x_U has to be Hermitian symmetric
if ( norm(x_U - flipud(conj(x_U))) > 1e-5 * norm(x_U) )
    error('x_U is not Hermitian symmetric!')
end

% Construct R_tilde, [1, Theorem 1]
R_tilde = toeplitz(x_U((length(x_U)+1)/2:end), x_U((length(x_U)+1)/2:-1:1));
% For numerical stability, compute the Hermitian part of R_tilde
R_tilde = (R_tilde + R_tilde')/2;

[EV, EW] = eig(R_tilde);
% Sort by the absolute values of eigenvalues
ew = diag(EW);
[~, II] = sort(abs(ew), 'descend');
% Sort eigenvectors
% U_all = [Us, Un]
% Us and Un are the signal and noise subspace, respectively.
U_all = EV(:, II);
% Check whether DD exceeds the dimension
if (2*DD+1 > length(x_U))
    disp('Warning, the number of sources exceeds the limit of x_U. Choose only one vector as noise subspace.');
    Un = U_all(:, end);
else
    Un = U_all(:, DD+1:end); % Noise subspace
end

% Dimension of the noise subspace
[D_1, Nullity] = size(Un);

if (nargout == 2 || nargout == 3)
    % Number of grid points
    Npt = 2^16;
    
    H = zeros(Npt, Nullity);
    [H(:, 1), w] = freqz(Un(:, 1), 1, Npt, 'whole');
    for ii = 2 : Nullity
        H(:, ii) = freqz(Un(:, ii), 1, Npt, 'whole');
    end
    theta_bar_grid = w / 2 / pi;
    theta_bar_grid(theta_bar_grid >= 1/2) = theta_bar_grid(theta_bar_grid >= 1/2) - 1;
    theta_bar_grid = fftshift(theta_bar_grid);
    
    % MUSIC spectrum
    P = zeros(Npt, 1);
    for kk = 1 : Npt
        P(kk) = 1/norm(H(kk, :))^2;
    end
    P = P / max(P);
    P = fftshift(P);
    varargout{1} = theta_bar_grid;
    varargout{2} = P;
end

if (nargout == 1 || nargout == 3)
    
    % Root MUSIC algorithm
    coef = zeros(1, 2*D_1-1);
    for ii = 1 : Nullity
        coef = coef + conv(Un(:, ii), conj(flipud(Un(:, ii)))).';
    end
    peak_locs = roots(coef);
    
    radius = abs(peak_locs);
    radius(radius > 1) = -radius(radius > 1);
    [~, II] = sort( radius, 'descend' );
    theta_bar_est = sort( angle(peak_locs(II(1:DD))).' / (2*pi) ).';
    
    if (nargout == 1)
        varargout{1} = theta_bar_est;
    else
        varargout{3} = theta_bar_est;
    end
    
end

if (nargout < 1 || nargout > 3)
    error('Output arguments are not supported!')
end

end

