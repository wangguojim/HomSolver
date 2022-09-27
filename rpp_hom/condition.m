function Con=condition(Q,pr,x,ld,ud)
% x=xcplex;ld=0;ud=1;
res=Q*x+pr;
eps=10^-10;
Jl=find(x<=ld+eps);
Jm=find(ld+eps<x&x<ud-eps);
Ju=find(x>=ud-eps);
Con=zeros(5,1);
if (~isempty(Jl))
    Con(1)=abs(min(min(x),0));
    Con(3)=abs(min(min(res(Jl)),0));
end
if (~isempty(Ju))
    Con(2)=max(max(x-ud),0);
    Con(5)=max(max(res(Ju)),0);
end
if (~isempty(Jm))
    Con(4)=max(abs(res(Jm)).*min(abs(x(Jm)),abs(1-x(Jm))));
end
