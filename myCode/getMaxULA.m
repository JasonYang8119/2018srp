function [ output_vec ] = getMaxULA( fieldArray )
%GETMAXULA ��ȡ��ľ���������ULA
%   ����ʵ�ʷ������У�ͨ�����õ��㷨�����ULA����
[n1,n2] = ndgrid(fieldArray);
n1_n2 = [-n1 - n2,n1 - n2,n1 + n2];
n1_n2_vec = unique(n1_n2(:), 'sorted');
% Locate the maximum contiguous set
D_MAX = max(n1_n2_vec);
D_full = (-D_MAX : D_MAX).';
H = setdiff(D_full, n1_n2_vec); % Holes
if (isempty(H))
    % D is hole-free
    output_vec = n1_n2_vec;
else
    N_max = min(abs(H)) - 1;
    output_vec = (-N_max : N_max).';
end
end

