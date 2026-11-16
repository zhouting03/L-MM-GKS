function S = lplq_cache(op,V,cap)
% 缓存未加权 V、AV、PsiV；FULL 模式后续按小批次增长。
k=size(V,2); alloc=min(cap,max(k,8));
S.V=zeros(op.n,alloc); S.Y=zeros(op.m,alloc); S.Z=zeros(op.r,alloc);
S.k=k; S.cap=cap; S.alloc=alloc; S.nCompress=0; S.nSafeguard=0;
S.maxPreserveError=0; S.maxDim=k; S.peakAllocatedMiB=8*(op.n+op.m+op.r)*alloc/2^20;
S.V(:,1:k)=V;
for j=1:k, S.Y(:,j)=op.A(V(:,j)); S.Z(:,j)=op.Psi(V(:,j)); end
end
