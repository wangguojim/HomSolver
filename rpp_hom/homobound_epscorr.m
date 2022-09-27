function [x,In_num,Out_num]=homobound_epscorr(Q,x,r,w,Jl,Jm,Ju,C,condQ) %% 算法存在问题收敛精度低
%该函数用来求解min_0=<x<=C 1/2x'Hx+r*x+w'x问题的解，其中min_x>=0 1/2x'Hx+r'*x的解已知


n=length(x);
epsilon=max(condQ/2^53,10^-8);   %设置为10^-8效果较好
tmax=1;
p=length(x);
t=0;
count=1;
trecord=0;
xrecord=x;
res=Q*sparse(x)+r;
In_num=0;
Out_num=0;

while(t<tmax)
    count        = count+1;
    fr=zeros(n,1);fr(Ju)=C;gr=Q*sparse(fr);
    if ~isempty(Jm)
    uv=Q(Jm,Jm)\(-[r(Jm)+gr(Jm),w(Jm)]);
    u=uv(:,1); v=uv(:,2);
    else
        u=[];v=[];
    end
    ux=zeros(n,1);
    vx=ux;ux(Jm)=u;vx(Jm)=v;
    uvg=Q*sparse([ux,vx]);
    ac=uvg(:,1)+r+gr;
    bd=uvg(:,2)+w;
    a=ac(Jl);b=bd(Jl);
    c=ac(Ju);d=bd(Ju);
    tall=-a./b;
    tupl=find(tall>trecord(end)+10*epsilon);
    [tjl,idl]=min(tall(tupl)) ;
    idl=tupl(idl);
    talu=-c./d;
    tupu=find(talu>trecord(end)+10*epsilon);
    [tju,idu]=min(talu(tupu)) ;
    idu=tupu(idu);
    tallx=-u./v;
    talux=(C-u)./v;
    tuplx=find(tallx>trecord(end)+10*epsilon);
    [tjlx,idlx]=min(tallx(tuplx)) ;
    idlx=tuplx(idlx);
    tupux=find(talux>trecord(end)+10*epsilon);
    [tjux,idux]=min(talux(tupux)) ;
    idux=tupux(idux);
    if(isempty(tjl))
        tjl=100;
    end
    if(isempty(tju))
        tju=101;
    end
    if(isempty(tjlx))
        tjlx=102;
    end
    if(isempty(tjux))
        tjux=103;
    end
    [tj,id]=min([tjl,tju,tjlx,tjux]);
    x=zeros(n,1);
    x(Jm)= u+min(tj,tmax)*v  ;
    res(Jl)=a+min(tj,tmax)*b;
    res(Ju)=c+min(tj,tmax)*d;
    
    if(min(x)<-100*epsilon)    
        id=0;
    end
    if(max(x)>C+100*epsilon)    
        id=0;
    end
    if(min(res(Jl))<-1000*epsilon)    
        id=0;
    end
    if(max(res(Ju))>1000*epsilon)      
        id=0;
    end
    %% ----------------------------------------------------------%
    if tj>=tmax
        tj=tmax;
    else
        if (id==1)
            Jlchu=Jl(idl);
            Jl=setdiff(Jl,Jlchu);
            Jm=union(Jm,Jlchu);
            In_num=In_num+1;
            tj=tjl;
        elseif (id==2)
            Juchu=Ju(idu);
            Ju=setdiff(Ju,Juchu);
            Jm=union(Jm,Juchu);
            In_num=In_num+1;
            tj=tju;
        elseif (id==3)
            Jljin=Jm(idlx);
            Jm=setdiff(Jm,Jljin);
            x(Jljin)=0;
            Jl=union(Jl,Jljin);
            Out_num=Out_num+1;
            tj=tjlx;
        elseif (id==4)
            Jujin=Jm(idux);
            Jm=setdiff(Jm,Jujin);
            x(Jujin)=C;
            Ju=union(Ju,Jujin);
            Out_num=Out_num+1;
            tj=tjux;
        elseif (id==0)
            f=r+tj*w;
            [Jl,Jm,Ju,exitflag]=eps_correct_bound(Q,f,Jl,Jm,Ju,C,epsilon);
            if exitflag==0  
                x=x.*(x>0&x<C)+C.*(x>=C);
                x(Jl)=0;x(Ju)=C;
                return;
            end
        end
    end
    
    x(Jl)=0;
    x(Ju)=C;
    xrecord(:,count)    = x;
    trecord(count) = tj ;
    t=tj ;
end