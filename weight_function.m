function [D, w] = weight_function( S, coarray_set )
%WEIGHT_FUNCTION: Returns the weight function of sensor array geometry S
%==============================================================================
%
%	Written by Chun-Lin Liu
%	E-mail: cl.liu@caltech.edu
%	Project website: 
%	http://systems.caltech.edu/dsp/students/clliu/coarray.html
%	Last revised on July 13, 2017
%
%   If you use this function, please cite related references [1-6].
%
%   The definition of weight function w(m) can be found in [1-6].
%   Input:
%       1) S: The discrete integer set
%       2) coarray_set: specify the output grid D. It has three options
%           'D': w is defined over D
%           'U': defined over U, the central maximum contiguous ULA portion
%           'V': defined over V, including the holes
%           The definitions of D, U, and V can be found in [6].
%
%   Output:
%       1) D: Output grid
%       2) w: weight function on D
%
%   Examples:
%       [D, w] = weight_function([0; 1; 4], 'D');
%       We obtain
%           D = [-4; -3; -1; 0; 1; 3; 4];
%           w = [ 1;  1;  1; 3; 1; 1; 1];
%
%       [D, w] = weight_function([0; 1; 4], 'U');
%       We obtain
%           D = [-1; 0; 1];
%           w = [ 1; 3; 1];
%
%       [D, w] = weight_function([0; 1; 4], 'V');
%       We obtain
%           D = [-4; -3; -2; -1; 0; 1; 2; 3; 4];
%           w = [ 1;  1;  0;  1; 3; 1; 0; 1; 1];
%
%	References
%	Journal
%   [1] C.-L. Liu and P. P. Vaidyanathan, â€œRemarks on the Spatial Smoothing Step in Coarray MUSIC,â€?IEEE Signal Processing Letters, vol. 22, no. 9, pp. 1438-1442, Sep. 2015.
%	[2] C.-L. Liu and P. P. Vaidyanathan, â€œSuper Nested Arrays: Linear Sparse Arrays with Reduced Mutual Coupling - Part I: Fundamentals,â€?IEEE Transactions on Signal Processing, vol. 64, no. 15, pp. 3997-4012, Aug. 2016. 
%	[3] C.-L. Liu and P. P. Vaidyanathan, â€œSuper Nested Arrays: Linear Sparse Arrays with Reduced Mutual Coupling - Part II: High-Order Extensions,â€?IEEE Transactions on Signal Processing, vol. 64, no. 16, pp. 4203-4217, Aug. 2016. 
%	Conference
%	[4] C.-L. Liu and P. P. Vaidyanathan, â€œSuper Nested Arrays: Sparse Arrays with Less Mutual Coupling than Nested Arrays,â€?in Proc. of 2016 IEEE International Conference on Acoustics Speech and Signal Processing (ICASSP 2016), pp. 2976-2980, Shanghai, China, Mar. 2016.
%	[5] C.-L. Liu and P. P. Vaidyanathan, â€œHigh Order Super Nested Arrays,â€?in Proc. of the Ninth IEEE Sensor Array and Multichannel Signal Processing Workshop (SAM 2016), Rio de Janeiro, Brazil, Jul. 2016.
%	[6] C.-L. Liu, P. P. Vaidyanathan and P. Pal, â€œCoprime Coarray Interpolation for DOA Estimation via Nuclear Norm Minimization,â€?in Proc. of 2016 IEEE International Symposium on Circuits and Systems (ISCAS 2016), pp. 2639-2642, Montreal, Canada, May 2016. 
% 
%==============================================================================

    % S is a vector
    if ( ~isvector(S) ||  ~isreal(S) )
        error('The parameter S has to be a a real vector')
    end
    
    S = unique(S);
    [n1, n2] = ndgrid(S);
    n1_n2_mat = n1 - n2;
    n1_n2_vec = n1_n2_mat(:);
    
    D_set = unique(n1_n2_vec, 'sorted');
    w_set = histc(n1_n2_vec, D_set);
    
    D_MAX = max(D_set); % The maximum element in the difference coarray
    
    
    switch(coarray_set)
        case 'D'
            D = D_set;
            w = w_set;
        case 'U'
            % Locate the maximum contiguous set
            D_full = (-D_MAX : D_MAX).';
            H = setdiff(D_full, D_set); % Holes
            if (isempty(H))
                % D is hole-free
                D = D_set; w = w_set;
            else
                N_max = min(abs(H)) - 1;
                D = (-N_max : N_max).';
                w = w_set(D_set >= -N_max & D_set <= N_max);
            end
        case 'V'
            % Locate the central maximum contiguous set
            D_full = (-D_MAX : D_MAX).';
            [~, ia, ~] = intersect(D_full, D_set);
            if (length(ia) == length(D_full))
                % D is hole-free
                D = D_set; w = w_set;
            else
                D = D_full;
                w = zeros(size(D));
                w(ia) = w_set;
            end
        case 'SD'
            %SD = negtive_sum + positive_sum +difference
            %difference have beed done: n1_n2_mat = n1 - n2
            %n1_n2_vec = n1_n2_mat(:) is a column vector
            %positive_sum
            n1_n2_psum = n1 + n2;
            %negtive_sum
            n1_n2_nsum = -n1 - n2;
            n1_n2_mat = [n1_n2_psum n1_n2_mat n1_n2_nsum];

            %new vector
            n1_n2_vec = n1_n2_mat(:);           
            D_MAX = max(abs(n1_n2_vec));
            D_full = (-D_MAX : D_MAX).';
            D_set = unique(n1_n2_vec, 'sorted');
            w_set = histc(n1_n2_vec, D_set);
            H = setdiff(D_full, D_set); % Holes
            if (isempty(H))
                % D is hole-free
                D = D_set; w = w_set;
            else
                N_max = min(abs(H)) - 1;
                D = (-N_max : N_max).';
                w = w_set(D_set >= -N_max & D_set <= N_max);
            end
        otherwise
            error('The parameter coarray_set is either ''D'', ''U'', ''SD'' ,or ''V''')
    end


end

