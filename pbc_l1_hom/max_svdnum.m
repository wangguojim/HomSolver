function max_svd=max_svdnum(A,chose)
err=1;
n=size(A,2);
x=randn(n,1);
norm_x=norm(x);
x=x/norm_x;
if chose==1
    while(err>1e-4)
        x1=A*x;
        norm_x1=norm(x1);
        err=abs(norm_x1-norm_x);
        norm_x=norm_x1;
        x=x1/norm_x1;
    end
    max_svd=norm_x1;
else
    while(err>1e-4)
        x1=A'*(A*x);
        norm_x1=norm(x1);
        err=abs(norm_x1-norm_x);
        norm_x=norm_x1;
        x=x1/norm_x1;
    end
    max_svd=sqrt(norm_x1);
end
