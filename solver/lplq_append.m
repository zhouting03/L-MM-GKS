function S = lplq_append(S,op,v)
% 只对新增方向做正向算子乘法；不重复计算已有列。
if S.k>=S.cap, error('LpLq:Capacity','Cannot append beyond planned capacity.'); end
if S.k==S.alloc
    newalloc=min(S.cap,S.alloc+8);
    S.V(:,newalloc)=0; S.Y(:,newalloc)=0; S.Z(:,newalloc)=0;
    S.alloc=newalloc;
end
S.k=S.k+1; S.V(:,S.k)=v;
S.Y(:,S.k)=op.A(v); S.Z(:,S.k)=op.Psi(v);
S.maxDim=max(S.maxDim,S.k);
S.peakAllocatedMiB=max(S.peakAllocatedMiB, ...
    8*(size(S.V,1)+size(S.Y,1)+size(S.Z,1))*S.alloc/2^20);
end
