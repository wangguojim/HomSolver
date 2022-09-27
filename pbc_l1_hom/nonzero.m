function nonz=nonzero(x)
n=length(x);
u=sort(abs(x),'descend');
normx=norm(x,1);
nonz=0;
for i=1:n
    if(u(i)==0)
        break;
    end
    sp=norm(u(1:i),1)/normx;
    if sp>0.999
        nonz=nonz+1;
        break;
    else
        nonz=nonz+1;
    end
end

