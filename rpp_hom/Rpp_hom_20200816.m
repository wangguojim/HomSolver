% % 该程序求解带盒式约束的非凸二次规划问题
% % min 1/2x'Qx+r'x
% %s.t.  0=<x<=ud
% % where Q is not positive definite


%% ----Random data------------------
n=1000;
global H Q;
A=randn(n,n);Q=A'+A+10*eye(n);issparse=0;
% A=sprandn(n,n,0.1);Q=A'+A+1*speye(n,n);issparse=1;
r=randn(n,1);

% %%
% r=pr;
n=length(r);
ld=0*ones(n,1);
ud=10*ones(n,1);  %upper bound

% 
delta=0.1;
tol=10^-11;tol2=2*10^-5;tol3=4*10^-5;
iterMax=100;
% % 
para.b=100;
%% --------------matlab、cplex、gurobi--------------------------%%
% tstart=tic;
% opt=optimoptions(@quadprog,'MaxIter',10000,'Display','iter','Algorithm','trust-region-reflective','TolFun',1e-40);
% tic;xmat=quadprog(Q,pr,[],[],[],[],zeros(n,1),ud*ones(n,1),[],opt);toc;  %CPLEXcplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub)
% telapsedmat = toc(tstart);
% fcn=@(x)(1/2*x'*Q*x+pr'*x);
%% ----------------------------------------------------------------%
tic;
DIS=0;
tstart=tic;
xk=zeros(n,1);
X=xk;
distance_err=1;
fval=[];normxk=0;
W=[];
res=Q*xk+r;
count=1;
while(distance_err>tol)
    
    Wl=find(xk<=ld&res<0);
    Wm=find(ld<xk&xk<ud);
    Wu=find(xk>=ud&res>0);
   
    if (distance_err/(normxk+1))>1e-2
        resWl=res(Wl);
        [~,ID1]=sort(resWl,'ascend');
        W1=Wl(ID1(1:min(length(Wl),para.b)));
        
        resWu=res(Wu);
        [~,ID2]=sort(resWu,'descend');
        W2=Wu(ID2(1:min(length(Wu),para.b)));
       
       
        
        W=union(Wm,union(W1,W2));
    else
        W=union(Wm,union(Wl,Wu));        
    end        
    lenW=length(W);
    QWW=Q(W,W);   
    rk=res(W)-QWW*xk(W);   
    EIG=(eig(full(QWW)));    
    ldW=ld(W);udW=ud(W);
    pos=delta-min(0,min(EIG));
    H=Q(W,W)+pos*speye(lenW);
    pr=rk-pos*xk(W);
    psk=xk(W);
    x0=zeros(lenW,1);
    vk=x0;
    t=1;
    iter=0;
    non_active_num=n;
    non_active_nump=non_active_num-1;
    Same_num=0;
    Same_Max=floor(100+n/200);
    err_xy=1;
    L=max(EIG)+pos;
    tol3=tol3/4;
    while(iter<iterMax && err_xy>tol3)
        tk=(1+sqrt(1+4*t^2))/2;
        gk=H*vk+pr;
        sk=(vk-1/L*gk);
        sk=ldW.*(sk<=ldW)+sk.*(sk>ldW&sk<udW)+udW.*(sk>=udW);
        vk=(sk+((t-1)/tk)*(sk-x0));
        t=tk;
        trun=tol2*norm(sk);
        non_active_num=sum(trun<sk&sk<W-trun);
        if(non_active_nump==non_active_num)
            Same_num=Same_num+1;
            non_active_nump=non_active_num;
        else
            Same_num=0;
            non_active_nump=non_active_num;
        end
        if(Same_num==Same_Max)
            break;
        end
        iter=iter+1;
        err_xy=norm(sk-x0)/norm(sk);
        x0=sk;
    end    
    
    if(norm(psk-sk)<10^-1)
        Jl=find(sk<tol);
        Jm=find(tol<sk&sk<udW-tol);
        Ju=find(sk>udW-tol);
        xinit=zeros(lenW,1);
        xinit(Jm)=sk(Jm);
        xinit(Ju)=min(ud);
        gk=H*xinit+pr;  %
        w=zeros(lenW,1);
        w(Jm)=-gk(Jm);
        w(Jl)=-min(gk(Jl))+10;
        w(Ju)=-max(gk(Ju))-10;      
%         sk1=parbdqp(xinit,pr+w,-w,Jl,Jm,Ju,min(ud));        
         [sk,~,~]=homobound_epscorr(H,xinit,pr+w,-w,Jl,Jm,Ju,min(ud),max(EIG)/delta);
      
    end
    normxk=sqrt(normxk^2+norm(sk)^2-norm(xk(W))^2);
    distance_err=norm(xk(W)-sk);
    if mod(count,100)==0
        res=Q*xk+r;
    else  
        hh=zeros(n,1);
        hh(W)=sk-xk(W);
        res=res+Q*sparse(hh);        
    end
    xk(W)=sk;
    count=count+1;
end

cc=[xk,res];
cc(W,:)
toc;
% telapsedAPG = toc(tstart);
% DISPPA=DIS(2:end);
% % figure
% % semilogy(1:length(DISPPA),DISPPA)
% %% ---------------------------------------------------------------
% % 
fcn=@(x)(1/2*x'*Q*x+r'*x);
grad=@(x)(Q*x+r);
fun     = @(x)fminunc_wrapper( x, fcn, grad);
opts    = struct( 'factr', 1e-22, 'pgtol', 1e-22, 'm', 20,'maxIts',4000);

opts.printEvery     = 100;
if n > 10000
    opts.m  = 50;
end
% Run the algorithm:
tstart=tic;[xlbfgsb, ~, info] = lbfgsb(fun, zeros(n,1), ud.*ones(n,1), opts );
telapsedlbfgsb = toc(tstart)
s1=fcn(xlbfgsb)
s2=fcn(xk)
APPA_convex
% 
% %% --------------------------------------------------------------%
% res=Q*xk+pr;
% Jm=find(10^-6<xk&xk<ud-10^-6);
% grad_APG=norm(res(Jm));
% res=Q*xmat+pr;
% Jm=find(10^-6<xmat&xmat<ud-10^-6);
% grad_TRF=norm(res(Jm));
% res=Q*xlbfgsb+pr;
% Jm=find(10^-6<xlbfgsb&xlbfgsb<ud-10^-6);
% grad_LBFGSB=norm(res(Jm));
% xAPG=xk;
% fprintf(1,'n=%d\n',n);
% fprintf(1,'PPA/APG-PAS time=%5.2f\t',telapsedAPG);
% fprintf(1,'grad_PPA=%5.2e\n',grad_APG);
% fprintf(1,'L-BFGSB time=%5.2f\t',telapsedlbfgsb);
% fprintf(1,'grad_LBFGSB=%5.2e\n',grad_LBFGSB);
% fprintf(1,'MATtime=%5.2f\t',telapsedmat);
% fprintf(1,'grad_TRF=%5.2e\n',grad_TRF);
% fprintf(1,'fPPA_LBFGSB=%5.2e\t',1/2*xAPG'*Q*xAPG+pr'*xAPG-(1/2*xlbfgsb'*Q*xlbfgsb+pr'*xlbfgsb));
% fprintf(1,'fPPA_TRF=%5.2e\t',1/2*xAPG'*Q*xAPG+pr'*xAPG-(1/2*xmat'*Q*xmat+pr'*xmat));
