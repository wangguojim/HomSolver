function C=randchoose(n,k,v)
% 随机选取k个间隔大于v的指标
% n=100
% v=4
% k=41
S=1:n;
C=[];
for i=1:k
    lenS=length(S);
    if lenS==0
        break
    end
    j=randperm(lenS,1);
    id=S(j);
    C(end+1)=id;
    S(max(1,j-v):min(j+v,lenS))=[];
end

    
    