function varargout = sample_covariance_to_difference_coarray( S, R_S_tilde, output_type )
%SAMPLE_COVARIANCE_TO_DIFFERENCE_COARRAY:
%   Convert the sample covariance matrix R_S_tilde to quantities on the
%   difference coarray.
%==============================================================================
%
%	Written by Chun-Lin Liu
%	E-mail: cl.liu@caltech.edu
%	Project website:
%	http://systems.caltech.edu/dsp/students/clliu/coarray.html
%	Last revised on July 13, 2017
%
%   If you use this function, please cite reference [1] or [2].
%
%   Input:
%
%       1) S: The discrete integer set for sensor locations
%       2) R_S_tilde: The sample covariance matrix defined on S
%       3) output_type: The quantities on the difference coarray [2].
%           This parameter is case-insensitive.
%           Please see [2] for the definitions of D, U, and V.
%
%           a) 'x_D': The sample autocorrelation vector on the difference
%           coarray D.
%           b) 'x_U': The sample autocorrelation vector on the ULA segment
%           of the difference coarray (set U).
%           c) 'x_V': The sample autocorrelation vector on the set V. Here
%           the coarray interpolation method is nuclear norm minimization.
%               Note: cvx package (http://cvxr.com/cvx/) is required here.
%           d) 'R_D': Augmented covariance matrix based on x_D.
%           e) 'R_U': Augmented covariance matrix based on x_U.
%           f) 'R_V': Augmented covariance matrix based on x_V.
%               Note: cvx package (http://cvxr.com/cvx/) is required here.
%
%   Output:
%
%       1) If there is only one output parameter, then it means either
%           x_D, x_U, x_V, R_D, R_U, or R_V.
%       2) If there are two output variables,
%           the first argument is either D, U, or V;
%           the second argument is either x_D, x_U, x_V, R_D, R_U, or R_V.
%
%   References
%   [1] C.-L. Liu and P. P. Vaidyanathan, 鈥淩emarks on the Spatial Smoothing Step in Coarray MUSIC,鈥?IEEE Signal Processing Letters, vol. 22, no. 9, pp. 1438-1442, Sep. 2015. 
%   [2] C.-L. Liu, P. P. Vaidyanathan and P. Pal, 鈥淐oprime Coarray Interpolation for DOA Estimation via Nuclear Norm Minimization,鈥?in Proc. of 2016 IEEE International Symposium on Circuits and Systems (ISCAS 2016), pp. 2639-2642, Montreal, Canada, May 2016. 
%
%==============================================================================

switch (upper(output_type))
    case 'X_D'
        
        [n1, n2] = ndgrid(S);
        [D, ~] = weight_function( S, 'D' );
        LEN_D = length(D);
        
        x_D = zeros( LEN_D, 1 );
        for mm = 1 : LEN_D
            % Select samples with the same lag
            data = R_S_tilde( n1 - n2 == D(mm) );
            % ====================
            % Average over all the samples with the same lag.
            % Definition 3 in [1]
            x_D(mm) = mean( data );
            % ====================
        end
        output_1 = D;
        output_2 = x_D;
        
    case 'X_U'
        
        [n1, n2] = ndgrid(S);
        [U, ~] = weight_function( S, 'U' );%得到最大长度的子序列
        LEN_U = length(U);
        
        x_U = zeros( LEN_U, 1 );
        for mm = 1 : LEN_U
            % Select samples with the same lag
            data = R_S_tilde( n1 - n2 == U(mm) );
            % ====================
            % Average over all the samples with the same lag.
            % Definition 3 in [1]
            x_U(mm) = mean( data );
            % ====================
        end
        output_1 = U;
        output_2 = x_U;
        
    case 'X_V'
        
        % Nuclear norm minimization
        [V, R_V] = sample_covariance_to_difference_coarray( S, R_S_tilde, 'R_V' );
        
        output_1 = V;
        output_2 = [R_V(1, end:-1:2).'; R_V(:, 1)];
        
    case 'R_D'
        
        [D, x_D] = sample_covariance_to_difference_coarray( S, R_S_tilde, 'x_D' );
        LEN_D = length(D);
        
        output_1 = D;
        output_2 = toeplitz(x_D((LEN_D+1)/2:end), x_D((LEN_D+1)/2:-1:1));
        
    case 'R_U'
        
        [U, x_U] = sample_covariance_to_difference_coarray( S, R_S_tilde, 'x_U' );
        LEN_U = length(U);
        
        output_1 = U;
        output_2 = toeplitz(x_U((LEN_U+1)/2:end), x_U((LEN_U+1)/2:-1:1));
        
    case 'R_V'
        
        % Sample covariance matrix on the difference coarray
        [D, x_D] = sample_covariance_to_difference_coarray( S, R_S_tilde, 'x_D' );
        
        % Nuclear norm minimization
        [V, ~] = weight_function( S, 'V' );
        LEN_V = length(V);
        Vp = V((1+LEN_V)/2:end);
        LEN_Vp = length(Vp);
        [~, index_S4p_pos, ~] = intersect( Vp, D(D >= 0) );
        cvx_begin % quiet
            variable R_V(LEN_Vp, LEN_Vp) hermitian toeplitz
            minimize( norm_nuc(R_V) )
            subject to
                R_V(index_S4p_pos, 1) == x_D(D >= 0);
        cvx_end
        
        output_1 = V;
        if ( strcmp(cvx_status, 'Infeasible') )
            output_2 = zeros(size(R_V));
        else
            output_2 = R_V;
        end
        
    otherwise
        error('output_type is invalid!')
end


if (nargout == 1)
    varargout{1} = output_2;
elseif (nargout == 2)
    varargout{1} = output_1;
    varargout{2} = output_2;
else
    error('Output arguments are not supported!')
end

end

