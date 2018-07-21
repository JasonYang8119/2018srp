function varargout = clscl_MUSIC( R_cov,S,DD,theta )
%COARRAY_MUSIC: MUSIC spectrum using sample cov
[EV, EW] = eig(R_cov);
% Sort by the absolute values of eigenvalues
ew = diag(EW);
[~, II] = sort(abs(ew), 'descend');
% Sort eigenvectors
% U_all = [Us, Un]
% Us and Un are the signal and noise subspace, respectively.
U_all = EV(:, II);
% Check whether DD exceeds the dimension
if (DD >= length(ew))
    disp('Warning, the number of sources exceeds the limit of x_U. Choose only one vector as noise subspace.');
    Un = U_all(:, end);
else
    Un = U_all(:, DD+1:end); % Noise subspace
end
%number of gird point
Npt = 1;
%search A_S
for ii=min(theta):0.1:max(theta)
    theta_sin = sin(pi*ii/180);
    a_theta = exp(1i * pi * S * theta_sin);
    P(Npt)=abs(1/((a_theta'*Un)*(a_theta'*Un)'));
    Npt = Npt + 1;
end
P = P/max(P);
varargout{1} = P;
varargout{2} = min(theta):0.1:max(theta);
end

