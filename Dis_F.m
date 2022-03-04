function [DATA_ORG] = Dis_F(data)
dist1 = pdist(data);
dist1=dist1';
k=1;

for i=1:size(data,1)
    for j=i+1:size(data,1)
        DATA_ORG(i,j)=dist1(k);
        k=k+1;
    end
end

DATA_ORG(i,:)=0;

for i=1:size(data,1)
    for j=i+1:size(data,1)
        DATA_ORG(j,i)=DATA_ORG(i,j);
    end
end

