function y = lplq_grad2(x,nx,ny,adjoint)
% 周期前向差分及其精确伴随，对应各向异性梯度正则。
if ~adjoint
    X=reshape(x,nx,ny);
    a=X([2:end,1],:)-X; b=X(:,[2:end,1])-X; y=[a(:);b(:)];
else
    n=nx*ny; a=reshape(x(1:n),nx,ny); b=reshape(x(n+1:2*n),nx,ny);
    X=a([end,1:end-1],:)-a+b(:,[end,1:end-1])-b; y=X(:);
end
end
