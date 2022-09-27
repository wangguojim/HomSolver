% % 该程序求解带盒式约束的非凸二次规划问题
% % min 1/2x'Qx+r'x
% %s.t.  0=<x<=ud
% % where Q is not positive definite
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----Random data------------------
% n=1000;
% % m=floor(0.1*n);
% global H Q;
% A=randn(n,n);Q=A'+A+10*eye(n,n);issparse=0;
% % A=sprandn(n,n,0.01);Q=A'+A+1*speye(n,n);issparse=1;
% pr=randn(n,1);
%% --------------real test problems-------------------



% ud=10;  %upper bound
% EIGH=eig(Q);
% enta=-min(EIGH);     %修正参数，使得问题变为凸问题
% delta=1;
tol=10^-11;tol2=4*10^-5;tol3=4*10^-4;
% iterMax=5*n;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L=max(EIGH)+enta+delta;
% % Q=A'*A-enta*eye(n,n);
%
% if issparse==1
%     H=Q+(enta+delta)*speye(n,n);
% else
%     H=Q+(enta+delta)*eye(n,n);
% end


%% --------------matlab、cplex、gurobi--------------------------%%



% tstart=tic;
% opt=optimoptions(@quadprog,'MaxIter',10000,'Display','iter','Algorithm','trust-region-reflective','TolFun',1e-150);
% tic;xmat=quadprog(Q,pr,[],[],[],[],zeros(n,1),ud*ones(n,1),[],opt);toc;  %CPLEXcplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub)
% telapsedmat = toc(tstart);

%% ----------------------------------------------------------------%
tstart=tic;
xk=zeros(n,1);
pxk=xk+1;
X=xk;
distance_err=1;
fval_x0=0;
DIS=[10,1];
OUTiter=0;
ac=0;
Act=0;
gamma0=0;
APPA_cho=1;
APPA_cho_max=1;  %1或者2,n越大，值越大
vk1=xk;Hxk=H*vk1;
count=1;
while(distance_err>tol)
    gamma=DIS(end)/DIS(end-1);
    if abs(gamma0-gamma)<1*10^-1&&gamma<1;
        APPA_cho=APPA_cho+1;
    else
        APPA_cho=0;
    end
    if APPA_cho==APPA_cho_max
        xk=xk+gamma/(1-gamma)*(xk-pxk);
                
        %num=600;xk=((1-gamma^num)*xk-(gamma-gamma^num)*pxk)/(1-gamma);
        xk=xk.*(xk>0&xk<ud)+ud*(xk>=ud);
        ac=0;
        APPA_cho=0;
    else
        ac=min(ac+1,APPA_cho_max);
    end
    r=pr-(delta+enta)*xk;
    pxk=xk;
    
    x0=xk;
    vk=x0;
    t=1;
    iter=0;
    non_active_num=n;
    non_active_nump=non_active_num-1;
    Same_num=0;
    Same_Max=floor(20+n/200);
    err_xy=1;
    while(iter<iterMax && err_xy>tol3)
        tk=(1+sqrt(1+4*t^2))/2;
         sk=vk-vk1;
         vk1=vk;
        if sum(sk~=0)<0.4*n;
            sk=sparse(sk);
        end
        count=count+1;
        if mod(count,3)==0 %2或3
            Hxk=H*vk;
        else
            Hxk=H*sk+Hxk;
        end
        
        gk=Hxk+r;
%          gk1=H*vk+r;        norm(gk-gk1)
        xk=(vk-1/L*gk);
        xk=xk.*(xk>0&xk<ud)+ud*(xk>ud);
        vk=(xk+((t-1)/tk)*(xk-x0));
        t=tk;
        trun=tol2*norm(xk);
        non_active_num=sum(trun<xk&xk<ud-trun);
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
        err_xy=norm(xk-x0)/norm(xk);
        x0=xk;
    end
%     fval_xk=1/2*xk'*Q*xk+pr'*xk;
    if(norm(pxk-xk)<10^-0)
        Jl=find(xk<tol);
        Jm=find(tol<xk&xk<ud-tol);
        Ju=find(xk>ud-tol);
        xinit=zeros(n,1);
        xinit(Jm)=xk(Jm);
        xinit(Ju)=ud;
        res=H*xinit+r;  %
        w=zeros(n,1);
        w(Jm)=-res(Jm);
        w(Jl)=-min(res(Jl))+10;
        w(Ju)=-max(res(Ju))-10;
        if issparse==1
            xk=parbdqpsparse(xinit,r+w,-w,Jl,Jm,Ju,ud);
        else
            xk=parbdqp(xinit,r+w,-w,Jl,Jm,Ju,ud);
        end
%         fval_xk=1/2*xk'*Q*xk+pr'*xk;
    end
   
    distance_err=norm(pxk-xk)
    DIS(end+1)=distance_err;
    X(:,end+1)=xk;
    OUTiter=OUTiter+1;
    Act(end+1)=length(Jm);
    gamma0=gamma;
end

telapsedAPG = toc(tstart);

DISAPPA=DIS(2:end);
figure;
semilogy(1:length(DISPPA),DISPPA,'k-.')
hold on;
semilogy(1:length(DISAPPA),DISAPPA,'k')
legend('PP-Hom','APP-Hom')
xlabel('Iter','FontSize',12)
ylabel('||x^{k}-x^{k+1}||','Interpreter','tex','FontSize',12)
% figure
% semilogy(1:length(DISAPPA),DISAPPA)
%% ---------------------------------------------------------------

fcn=@(x)(1/2*x'*Q*x+pr'*x);
grad=@(x)(Q*x+pr);
fun     = @(x)fminunc_wrapper( x, fcn, grad);
opts    = struct( 'factr', 1e-191, 'pgtol', 1e-1120, 'm', 50,'maxIts',10000);

opts.printEvery     = 100;
if n >= 10000
    opts.m  = 50;
end
%%Run the algorithm:
tstart=tic;[xlbfgsb, ~, info] = lbfgsb(fun, zeros(n,1), ud*ones(n,1), opts );
telapsedlbfgsb = toc(tstart)

%% --------------------------------------------------------------%
res=Q*xk+pr;
Jm=find(10^-6<xk&xk<ud-10^-6);
grad_APG=norm(res(Jm));
res=Q*xmat+pr;
Jm=find(10^-6<xmat&xmat<ud-10^-6);
grad_TRF=norm(res(Jm));
res=Q*xlbfgsb+pr;
Jm=find(10^-6<xlbfgsb&xlbfgsb<ud-10^-6);
grad_LBFGSB=norm(res(Jm));
xAPG=xk;
fprintf(1,'n=%d\n',n);
fprintf(1,'APPA/APG-PAS time=%5.2f\t',telapsedAPG);
fprintf(1,'grad_APPA=%5.2e\n',grad_APG);
fprintf(1,'L-BFGSB time=%5.2f\t',telapsedlbfgsb);
fprintf(1,'grad_LBFGSB=%5.2e\n',grad_LBFGSB);
fprintf(1,'MATtime=%5.2f\t',telapsedmat);
fprintf(1,'grad_TRF=%5.2e\n',grad_TRF);
% fprintf(1,'fAPPA_LBFGSB=%5.2e\t',1/2*xAPG'*Q*xAPG+pr'*xAPG-(1/2*xlbfgsb'*Q*xlbfgsb+pr'*xlbfgsb));
% fprintf(1,'fAPPA_TRF=%5.2e\t',1/2*xAPG'*Q*xAPG+pr'*xAPG-(1/2*xmat'*Q*xmat+pr'*xmat));

