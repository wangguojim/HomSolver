% % 该程序求解带盒式约束的非凸二次规划问题
% % min 1/2x'Qx+r'x
% %s.t.  0=<x<=ud
% % where Q is not positive definite


%% ----Random data------------------

% n=10000;
% 
%  %
% issparse=0;
% % 
% % 
% if issparse==1
%             A=sprandn(n,n,sparsity);Q=A'+A+1*speye(n,n);%10万.0.01~0.001，统计Q稀疏度：sum(sum(Q~=0))/n/n
%         else
%             A=randn(n,n);Q=A'+A+10*eye(n); %2,3万
%         end
%         r=randn(n,1);pr=r;
% clear A;
% r=randn(n,1);
% pr=r;
% n=length(r);
fcn=@(x)(1/2*x'*Q*x+r'*x);

ud=10;
lb=0*ones(n,1);
ub=10*ones(n,1);  %upper bound

delta=0.01;
tol=10^-11;tol2=2*10^-5;tol3=4*10^-5;
iterMax=100;
% %
if n<100000
    para.b=100;
    para.b1=1000;
else
    para.b=500;
    para.b1=5000;
end
%% ----------------------------------------------------------------%
time=tic;
DIS=0;
tstart=tic;
xk=zeros(n,1);
X=xk;
distance_err=1;
fval=[];normxk=0;
W=[];
res=Q*xk+r;
count=1;
EigRec=[];
res_err=1;
DISARPP=[10,1];
APPA_cho_max=2;  %1或者2,n越大，值越大
DIS=[10,1];
OUTiter=0;
ac=0;
Act=0;
gamma0=0;
POS=[];

while(res_err>tol )
    gamma=DISARPP(end)/DISARPP(end-1);
    if abs(gamma0-gamma)<1*10^-1&&gamma<1;
        APPA_cho=APPA_cho+1;
    else
        APPA_cho=0;
    end
    if APPA_cho==APPA_cho_max
        
        
        
%         
        xW1=xk(W);
        xW2=xW1+gamma/(1-gamma)*(xW1-pxW);
        
        xW2= xW2.*( xW2>0& xW2<ud)+ud*( xW2>=ud);
        xW2_W1=xW2-xW1;
        hs=zeros(n,1);
        hs(W)=xW2_W1;
        res=res+Q*sparse(hs);
        x1=xk;
        fvalue1=fcn(x1);
        x2=xk;
        
        x2(W)=xW2;
        fvalue2=fcn(x2);
        fvalue_err=fvalue2-fvalue1;
        if fvalue_err>0
            
            xk(W)=xW2;

        end
%         norm(res-Q*xk-r)
        %         xk1_vk=xk-vk1;
        %         if sum(xk1_vk~=0)<0.4*n;
        %             xk1_vk=sparse(xk1_vk);
        %         end
        %         Qxk1=Q*xk1_vk+Hvk-(enta+delta)*vk1;
        %
        %
        %         fvaule=(0.5*Qxk1+pr)'*xk;
        %
        %         %num=600;xk=((1-gamma^num)*xk-(gamma-gamma^num)*pxk)/(1-gamma);
        %         xk2=xk2.*(xk2>0&xk2<ub)+ub*(xk2>=ub);
        %         y=xk2-xk1;
        %
        %         if sum(y~=0)<0.4*n;
        %             y=sparse(y);
        %         end
        %         Qy=Q*y;
        %
        %
        %         coex=y'*Qy;
        %         coey=(Qxk1+pr)'*y;
        %         minp=-coey/coex;
        %         if coex>0
        %             cut=max(0,min(1,2*minp));
        %         else
        %             [~,aj]=min([abs(1-minp),abs(minp)]);
        %             if aj==1
        %                 cut=1;
        %
        %             else
        %                 cut=0;
        %             end
        %         end
        %
        %          fvdesc=coex/2+coey;
        % %          norm(fvdesc-(fx(xk2)-fx(xk1)))
        %         if fvdesc<0
        %             xk=xk2;
        %         else
        %             xk=xk1+cut*y;
        %         end
        %    xk=xk2;
        
        ac=0;
        APPA_cho=0;
    else
        ac=min(ac+1,APPA_cho_max);
    end
    Wl=find(xk<=lb&res<0);
    Wm=find(lb<xk&xk<ub);
    Wu=find(xk>=ub&res>0);
    
    if (distance_err/(normxk+1))>1e-3 || length(Wl)>para.b1/10  || length(W2)>para.b1/10 || length(Wm)>(para.b1)
        resWl=res(Wl);
        [~,ID1]=sort(resWl,'ascend');
        W1=Wl(ID1(1:min(length(Wl),para.b)));
        
        resWu=res(Wu);
        [~,ID2]=sort(resWu,'descend');
        W2=Wu(ID2(1:min(length(Wu),para.b)));
        
        resWm=res(Wm);
        [~,ID3]=sort(abs(resWm),'descend');
        
        W3=Wm(ID3(1:min(para.b1,length(Wm))));
        W=union(W3,union(W1,W2));
        
    else
        W1=Wl;W2=Wu;
        W=union(Wm,union(Wl,Wu));
    end
    
    res_err1=0;  res_err2=0;  res_err3=0;
    if ~isempty (Wm)
        res_err1=max(abs(res(Wm)));
    end
    if ~isempty (W1)
        res_err2=min(res(W1));
    end
    if ~isempty (W2)
        res_err3=max(res(W2));
        
    end
    res_err=abs(res_err1)+abs(res_err2)+abs(res_err3)
     
    if length(W)==0
        break;
    end
    lenW=length(W);
    QWW=Q(W,W);
    rk=res(W)-QWW*xk(W);
    %EIG=(eig(full(QWW)));
    
    [mineig,maxabs_eig]=mineigs(QWW,1e-4);
    %    EigRec(end+1)= mineig-0.01;
    
    ldW=lb(W);udW=ub(W);
    pos=delta-min(0,mineig-0.01);
    POS(end+1)=delta-pos;
    H=QWW+pos*speye(lenW);
    pk=rk-pos*xk(W);
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
    Lk=maxabs_eig+pos;
    tol3=tol3/4;
    while(iter<iterMax && err_xy>tol3)
        tk=(1+sqrt(1+4*t^2))/2;
        gk=H*vk+pk;
        sk=(vk-1/Lk*gk);
        sk=ldW.*(sk<=ldW)+sk.*(sk>ldW&sk<udW)+udW.*(sk>=udW);
        vk=(sk+((t-1)/tk)*(sk-x0));
        t=tk;
        trun=tol2*norm(sk);
        non_active_num=sum(trun<sk&sk<ub(W)-trun);
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
    
    if(norm(psk-sk)<10^1)
        Jl=find(sk<tol);
        Jm=find(tol<sk&sk<udW-tol);
        Ju=find(sk>udW-tol);
        xinit=zeros(lenW,1);
        xinit(Jm)=sk(Jm);
        xinit(Ju)=min(ub);
        gk=H*xinit+pk;  %
        w=zeros(lenW,1);
        w(Jm)=-gk(Jm);
        w(Jl)=-min(gk(Jl))+10;
        w(Ju)=-max(gk(Ju))-10;
        %         sk1=parbdqp(xinit,r+w,-w,Jl,Jm,Ju,min(ud));
        [sk,~,~]=homobound_epscorr(H,xinit,pk+w,-w,Jl,Jm,Ju,min(ub),L/delta);
        
    end
    % fvaue_decresase=1/2*sk'*QWW*sk+pk'*sk-1/2*xk(W)'*QWW*xk(W)+pk'*xk(W)
    normxk=sqrt(normxk^2+norm(sk)^2-norm(xk(W))^2);
    distance_err=norm(xk(W)-sk);
    pxW=xk(W);
    if mod(count,100)==0
        xk(W)=sk;
        res=Q*xk+r;
    else
        hh=zeros(n,1);
        hh(W)=sk-xk(W);
        res=res+Q*sparse(hh);
        xk(W)=sk;
    end
    count=count+1;
    gamma0=gamma;
    DISARPP(end+1)=distance_err;
end
time=toc(time);
cc=[xk,res];


fcn=@(x)(1/2*x'*Q*x+pr'*x);
fvalue=fcn(xk)
res=Q*xk+r;
Jm=find(lb+1e-6<xk&xk<ub-1e-6);
grad_APG=norm(res(Jm))
DISARPP=DISARPP(2:end);
% figure
semilogy(1:length(DISARPP),DISARPP)
% %% ---------------------------------------------------------------
% %
% fcn=@(x)(1/2*x'*Q*x+r'*x);
% grad=@(x)(Q*x+r);
% fun     = @(x)fminunc_wrapper( x, fcn, grad);
% opts    = struct( 'factr', 1e-22, 'pgtol', 1e-22, 'm', 20,'maxIts',4000);
%
% opts.printEvery     = 100;
% if n > 10000
%     opts.m  = 50;
% end
% % Run the algorithm:
% tstart=tic;[xlbfgsb, ~, info] = lbfgsb(fun, zeros(n,1), ud*ones(n,1), opts );
% telapsedlbfgsb = toc(tstart)
%

%% --------------matlab、cplex、gurobi---------trr-----------------%%
% tstart=tic;
% opt=optimoptions(@quadprog,'MaxIter',10000,'Display','iter','Algorithm','trust-region-reflective','TolFun',1e-40);
% tic;xmat=quadprog(Q,r,[],[],[],[],zeros(n,1),ud,[],opt);toc;  %CPLEXcplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub)
% telapsedmat = toc(tstart);
%
% fvalue=fcn(xmat)
% res=Q*xmat+r;
% Jm=find(10^-6<xmat&xmat<ud-10^-6);
% grad_TRF=norm(res(Jm))
%% LBFGS
% res=Q*xlbfgsb+r;
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
% fprintf(1,'fPPA_LBFGSB=%5.2e\t',1/2*xAPG'*Q*xAPG+r'*xAPG-(1/2*xlbfgsb'*Q*xlbfgsb+r'*xlbfgsb));
% fprintf(1,'fPPA_TRF=%5.2e\t',1/2*xAPG'*Q*xAPG+r'*xAPG-(1/2*xmat'*Q*xmat+r'*xmat));
