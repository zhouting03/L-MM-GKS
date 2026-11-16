function V = lplq_keep_basis(mandatory,history,maxcols,tol)
% 先保留当前解/梯度；归一后双重 MGS，零列与相关列不占位置。
n=size(mandatory,1);
if isempty(mandatory), n=size(history,1); end
V=zeros(n,maxcols); k=0;
for group=1:2
    if group==1, X=mandatory; else, X=history; end
    for j=1:size(X,2)
        v=X(:,j); nv=norm(v);
        if nv==0, continue; end
        v=v/nv;
        if k>0
            Q=V(:,1:k); v=v-Q*(Q'*v); v=v-Q*(Q'*v);
        end
        nv=norm(v);
        if nv>tol
            k=k+1; V(:,k)=v/nv;
            if k==maxcols, V=V(:,1:k); return; end
        end
    end
end
V=V(:,1:k);
end
