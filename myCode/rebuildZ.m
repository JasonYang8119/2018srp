function outputZ = rebuildZ( R_all,vecULA,fieldArray )
%REBUILDZ 从三个处理后的矩阵中rebuildZ
%   R1是相减，R2正相加，R3负相加，输出形如vecULA
outputZ = zeros(1,length(vecULA));
[n1,n2] = ndgrid(fieldArray);
n1_n2 = [n1 - n2,n1 + n2,-n1 - n2];

for i = 1:length(vecULA)
    tempData = R_all(n1_n2 == vecULA(i));
    outputZ(i) = mean(tempData);
end
end

