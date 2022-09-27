function x=parbdqp_real(x,r,w,Jl,Jm,Ju,ud) %% 算法存在问题收敛精度低
%该函数用来求解min_x>=0 1/2x'Hx+r*x+w'x问题的解，其中min_x>=0 1/2x'Hx+r'*x的解已知
if(isempty(Jm))
    return;
end
%% -------------------------------
global H;
n=length(x);
eps=10^-8;   %设置为10^-8效果较好
tmax=1;
p=length(x);
t=0;
count=1;
trecord=0;
xrecord=x;
while(t<tmax)  
    count        = count+1;      
%     u            = -solvep(R,r(Jm)+ud*sum(H(Jm,Ju),2)) ;
%     v            = -solvep(R,w(Jm)) ;
    fr=zeros(n,1);
    fr(Ju)=ud;
    gr=H*sparse(fr);
%     uv=-solvep(R,[r(Jm)+gr(Jm),w(Jm)]); 
    uv=H(Jm,Jm)\-[r(Jm)+gr(Jm),w(Jm)];
    u=uv(:,1);
    v=uv(:,2); 
    tj      = tmax   ;
    jinclude           =0;   
    jinid=-1;
    chuid=-1;
     % find the additional index 
    ux=zeros(n,1);    
    vx=ux;    
    ux(Jm)=u;
    vx(Jm)=v;
    uvg=H*sparse([ux,vx]); 
    ac=uvg(:,1)+r+gr;
    bd=uvg(:,2)+w;
    a=ac(Jl);
    b=bd(Jl);
    c=ac(Ju);
    d=bd(Ju);    
    talld=-a./b;
    tupid=find(talld>trecord(end)+eps);
    [tjl,idl]=min(talld(tupid)) ;
    talud=-c./d;
    tupud=find(talud>trecord(end)+eps);
    [tju,idu]=min(talud(tupud)) ;
    if(isempty(tjl))
        tjl=100;
    end
    if(isempty(tju))
        tju=101;
    end   
    if(tjl>=1)&&(tju>=1)
        tj=tmax;
    else
        if (tjl<tju)
            jinclude=Jl(tupid(idl));
            tj=tjl;
            jinid=0;
        else
            jinclude=Ju(tupud(idu));
            tj=tju;
            jinid=1;
        end
    end
    x(Jm)          = u+tj*v  ;
    %% -------------------------------------------------------%%
    % find the deletion index
    lzero=x(Jm);
    LJm=length(lzero);
    idm=0;
    tjc=tj;   
    for l=1:LJm  
        if (lzero(l)<-eps)            
            tb=trecord(end)+(tjc-trecord(end))*xrecord(Jm(l),end)/abs(xrecord(Jm(l),end)-x(Jm(l)));            
            if tb<tj&&tb>(trecord(end)+eps)
                tj=tb;
                idm=l;
                jinid=-1;
                chuid=0;
            end
        elseif  (lzero(l)>ud+eps)
            tb=trecord(end)+(tjc-trecord(end))*abs(xrecord(Jm(l),end)-ud)/abs(xrecord(Jm(l),end)-x(Jm(l)));            
            if tb<tj&&tb>(trecord(end)+eps)
                tj=tb;
                idm=l;
                jinid=-1;
                chuid=1;
            end            
        end
    end
    
    if(abs(tj-tmax)<eps)
        tj=tmax;
    end 
     x(Jm)          = u+tj*v;
     
    %% -----------------------------------------------------------%
   % flash the active set
    if(jinid==-1)
        if(chuid==0)
            Jl(end+1)=Jm(idm);
            x(Jm(idm))=0;
            Jm=setdiff(Jm,Jm(idm));
            
        elseif(chuid==1)
            Ju(end+1)=Jm(idm);
            x(Jm(idm))=ud;
            Jm=setdiff(Jm,Jm(idm));
            
        end
    elseif(jinid==0)
        Jm(end+1)=jinclude;
        Jl=setdiff(Jl,jinclude);
    else
        Jm(end+1)=jinclude;
        Ju=setdiff(Ju,jinclude);
    end
        
    %% ----------------------------------------------------------%     
    xrecord(:,count)    = x;    
    trecord(count) = tj ;  
    t=tj ;
end