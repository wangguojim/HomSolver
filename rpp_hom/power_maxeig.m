function [maxabs_eig,uv]=power_maxeig(M,uu,tol)
err=1;
Hc1=10^9;
while(err>tol)
    uv=M*uu;
    Hc=(norm(uv));
    uv=uv/Hc;
    err=abs(Hc-Hc1);
    Hc1=Hc;
    uu=uv;
end
uv=M*uu;
[~,j]=max(abs(uu));
if abs(uu(j))<1e-12
    maxabs_eig=0;
else
    maxabs_eig=uv(j)/uu(j);
end
end