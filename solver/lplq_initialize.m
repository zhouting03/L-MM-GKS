function [S,st,used] = lplq_initialize(op,d,o,st,cap)
% 两阶段初始化与论文式流程对应，首次代理求解计入总预算。
used=0;
if ~isempty(o.V0)
    V=lplq_keep_basis([st.x,st.g],o.V0,o.kmin,o.orthTol);
    S=lplq_cache(op,V,cap); return;
end
ell=min(o.initDim,o.kmax-2);
seed=lplq_gkb_seed(op,d,st.wd,ell,o.orthTol);
V=lplq_keep_basis([st.x,st.g],seed,min(cap,ell+2),o.orthTol);
S=lplq_cache(op,V,cap);
if strcmp(o.init,'simple'), return; end
[st,meta]=lplq_step(S,st,op,d,o); used=1;
initPeak=S.maxDim; initAlloc=S.peakAllocatedMiB;
seed=lplq_krylov_seed(op,st,d,o,o.kmin-2);
V=lplq_keep_basis([st.x,st.g],seed,o.kmin,o.orthTol);
S=lplq_cache(op,V,cap);
S.maxDim=max(S.maxDim,initPeak); S.peakAllocatedMiB=max(S.peakAllocatedMiB,initAlloc);
S.nSafeguard=double(meta.safeguard);
end
