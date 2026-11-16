function S = lplq_compress(S,st,op,d,o)
% 原文式保留 [VW,x,g]：优先保留 x,g，当前权重下重建当前列数的小矩阵。
oldk=S.k; V=S.V(:,1:oldk);
if strcmp(o.mode,'RESTART')
    keep=2; W=zeros(oldk,0);
else
    keep=o.kmin;
    [~,fac]=lplq_projected_qr(S,st,d,o);
    a=V'*st.x;
    W=lplq_select_basis(fac,a,min(keep-2,oldk),o,S.nCompress+1);
end
newV=lplq_keep_basis([st.x,st.g],V*W,keep,o.orthTol);
T=V'*newV; defect=norm(newV-V*T,'fro')/max(1,norm(newV,'fro'));
if defect<=o.cacheTol
    newY=S.Y(:,1:oldk)*T; newZ=S.Z(:,1:oldk)*T;
else
    % 若梯度不完全在旧空间，直接调用算子，避免投影掉空间外分量。
    newY=zeros(op.m,size(newV,2)); newZ=zeros(op.r,size(newV,2));
    for j=1:size(newV,2), newY(:,j)=op.A(newV(:,j)); newZ(:,j)=op.Psi(newV(:,j)); end
end
k=size(newV,2); S.V(:,1:k)=newV; S.Y(:,1:k)=newY; S.Z(:,1:k)=newZ;
S.V(:,k+1:oldk)=0; S.Y(:,k+1:oldk)=0; S.Z(:,k+1:oldk)=0;
S.k=k; S.nCompress=S.nCompress+1;
ex=norm(st.x-newV*(newV'*st.x))/max(norm(st.x),realmin);
eg=norm(st.g-newV*(newV'*st.g))/max(norm(st.g),realmin);
S.maxPreserveError=max(S.maxPreserveError,max(ex,eg));
if max(ex,eg)>1e-8, error('LpLq:Preserve','Compression lost the solution or complete gradient.'); end
end
