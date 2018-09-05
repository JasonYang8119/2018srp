function outputZ = rebuildZ( R_all,vecULA,fieldArray )
%REBUILDZ �����������ľ�����rebuildZ
%   R1�������R2����ӣ�R3����ӣ��������vecULA
outputZ = zeros(1,length(vecULA));
[n1,n2] = ndgrid(fieldArray);
n1_n2 = [n1 - n2,n1 + n2,-n1 - n2];

for i = 1:length(vecULA)
    tempData = R_all(n1_n2 == vecULA(i));
    outputZ(i) = mean(tempData);
end
end

