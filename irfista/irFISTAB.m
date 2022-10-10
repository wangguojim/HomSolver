% Reduced FISTA

%% ----------------------testing-----------------------

% ----------------SVC data set
dataId =5;
switch dataId
    case 1,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\covtype.libsvm.binary.scale.mat')
    case 2,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\avazu-app.mat')
    case 3,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\epsilon_normalized.t.mat')
    case 4,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\HIGGS.mat')
    case 5,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\kdda.mat')
    case 6,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\kddb.mat')
    case 7,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\news20.binary.mat')
    case 8,
        
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\rcv1_test.binary.mat')
    case 9,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\url_combined.mat')
    case 10,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\webspam.mat')
    case 11,
        train_data = load('F:\wang\study\program\program\lasso\pbcm  homotopy for lasso\data\mnist.mat')
        
        
        
end
% %%
A=train_data.data{1}(:,1:min(100000000,length(train_data.data{2}))-17);Y=train_data.data{2}(1:size(A,2));
sample_id=length(train_data.data{2})-1;y=full(train_data.data{1}(:,sample_id));%+1*randn(size(A,1),1);
% -------------------random dataset-------------------%
% m=1000;n=10000;
% A=randn(m,n);y=randn(m,1);
%%

% % %%
[m,n]=size(A);
LambdaMax=max(abs(A'*y));
para.lamb=1e-3*LambdaMax;
kt=@(x) (A'*(A*x-y));kt1=@(x)(abs(kt(x))-para.lamb);kkt=@(x) (x~=0).*(kt(x)+para.lamb*sign(x))+kt1(x).*(kt1(x)>0).*(x==0);

%%

x=zeros(n,1);
err=1;
iter=1;
para.itermax=100;
para.ITERMAX=1000;
para.OutIterMax=300;
para.errtol=1e-1;
para.KKTErrtol=1e-7;
para.t=1;
para.zero_add=min(m,3000);
kktcond=1;
para.chose=2;
% L=max_svdnum(A,2)^2;%
% x=FISTA(A,y,para.lamb,L,1e-10,20);
%%  µÍ¾«¶È½â
time=tic;
Ax=A*x;
time1=[]
fvalue1=[];
len_W=[];
Liter=[];
Nzero=[]
while(kktcond>para.KKTErrtol&&iter<para.OutIterMax)
    gk=A'*(Ax-y);
    reswhole=gk+para.lamb*sign(x);
    W0=find(x==0&abs(gk)>(0.9999*para.lamb));
    [~,W0sort]=sort(abs(gk(W0)));
    W1=W0(W0sort(end-min(length(W0),para.zero_add)+1:end));
    W2=find(x~=0);
    W=union(W1,W2);
    kktcond=norm(reswhole(W2))+norm(abs(gk(x==0&abs(gk)>(para.lamb)))-para.lamb)
    AW=A(:,W);
    para.L=1*max_svdnum(AW,2)^2;
    para.itermax=min(para.ITERMAX,para.itermax*1.5);
    para.errtol=max(min(para.errtol/10,kktcond/10),para.KKTErrtol^2);%para.errtol=min(para.errtol/10,1e-11);
       
    if para.chose==1
        u =subFISTA(AW,y,x(W),para);
    else
        HWW=full(AW'*AW);r=AW'*(-y);
        u =subFISTA(HWW,r,x(W),para);
    end
    time1(end+1)=toc(time);
    fvalue1(end+1)=0.5*norm(AW*u-y)^2+para.lamb*sum(abs(u));
    Ax=Ax+AW*(u-x(W));
    x(W)=u;
    iter=iter+1;
    len_W(end+1)=length(W);
    Liter(end+1)=para.L;
    Nzero(end+1)=length(W2);
end
time=toc(time);
xrf=x;
kkt_xrf=norm(kkt(xrf))
f=@(x) (1/2*norm(A*x-y)^2+para.lamb*sum(abs(x)));
 
