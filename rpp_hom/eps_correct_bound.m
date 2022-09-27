function  [Jl,Jm,Ju,exitflag]=eps_correct_bound(Q,f,Jl,Jm,Ju,C,epsilon)
%epslion correction to make the solution satisfy the KKT conditions in e-precision

n=length(f);
restol=1000*epsilon;xtol=100*epsilon;
loop=1;
x=zeros(n,1); x(Jl)=0;x(Ju)=C; 
iter=0;
exitflag=1;
while(loop) 
    iter=iter+1;
    if iter>n
        exitflag=0;
        return;
    end
    fr=zeros(n,1);fr1=fr;
    fr(Ju)=C;gr=Q*sparse(fr);    
    x(Jm)=Q(Jm,Jm)\(-f(Jm)-gr(Jm));
    fr1(Jm)=x(Jm);
    res=Q*sparse(fr1)+gr+f;
    if min(x)<min(-xtol,-1e-5)
        [~,idml]=min(x(Jm));
        Jljin=Jm(idml);
        Jm=setdiff(Jm,Jljin);
        x(Jljin)=0;
        Jl=union(Jl,Jljin);
        continue;
    end
    if max(x)>(C+max(xtol,1e-5))
        [~,idmu]=max(x(Jm));
        Jujin=Jm(idmu);
        Jm=setdiff(Jm,Jujin);
        x(Jujin)=C;
        Ju=union(Ju,Jujin);
        continue;
    end
    if min(res(Jl))<min(-restol,-1e-5)
        [~,idlm]=min(res(Jl));
        Jlchu=Jl(idlm);
        Jl=setdiff(Jl,Jlchu);
        Jm=union(Jm,Jlchu);
        continue;
    end
    if max(res(Ju))>max(restol,1e-5)
        [~,idum]=max(res(Ju));
        Juchu=Ju(idum);
        Ju=setdiff(Ju,Juchu);
        Jm=union(Jm,Juchu);
        continue;
    end
    loop=0;
end

