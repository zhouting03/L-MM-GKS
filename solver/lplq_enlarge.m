function [S,st,records] = lplq_enlarge(S,st,op,d,o,budget,target,startIndex)
% 每个循环解一次代理、更新双权重，再扩展最新完整梯度。
records=cell(budget,1); used=0;
for j=1:budget
    [next,meta]=lplq_step(S,st,op,d,o); st=next;
    if mod(startIndex+j,o.refreshEvery)==0
        st=lplq_state(op,d,st.x,o); % 用真实算子刷新缓存组合误差。
    end
    if meta.safeguard
        old=S; V=lplq_keep_basis([st.x,st.g],old.V(:,1:old.k),o.kmin,o.orthTol);
        S=lplq_cache(op,V,old.cap); S.nCompress=old.nCompress;
        S.nSafeguard=old.nSafeguard+1; S.maxDim=max(old.maxDim,S.k);
        S.peakAllocatedMiB=max(old.peakAllocatedMiB,S.peakAllocatedMiB);
        S.maxPreserveError=old.maxPreserveError;
    elseif st.gnorm>target
        Q=S.V(:,1:S.k); v=st.g-Q*(Q'*st.g); v=v-Q*(Q'*v); nv=norm(v);
        if nv>o.orthTol*max(st.gnorm,realmin) && S.k<S.cap
            S=lplq_append(S,op,v/nv);
        end
    end
    used=used+1; r=lplq_record(S,st,o,op);
    r.relStep=meta.relStep; r.rankDeficient=double(meta.rankDeficient); records{used}=r;
    if st.gnorm<=target, break; end
    if ~strcmp(o.mode,'FULL') && S.k>=o.kmax, break; end
end
records=records(1:used);
end
