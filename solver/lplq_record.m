function r = lplq_record(S,st,o,op)
% 只记录标量历史；不存每一步的高维解或基。
r.objective=st.F; r.gradNorm=st.gnorm; r.rre=NaN;
if ~isempty(o.xTrue), r.rre=norm(st.x-o.xTrue)/max(norm(o.xTrue),realmin); end
r.basisDim=S.k; r.peakBasisDim=S.maxDim;
r.activeCacheMiB=8*(size(S.V,1)+size(S.Y,1)+size(S.Z,1))*S.k/2^20;
r.allocatedCacheMiB=8*(numel(S.V)+numel(S.Y)+numel(S.Z))/2^20;
r.peakAllocatedCacheMiB=S.peakAllocatedMiB;
r.nCompress=S.nCompress; r.nSafeguard=S.nSafeguard;
r.preserveError=S.maxPreserveError;
c=op.getCounts(); r.callsA=c.A; r.callsAt=c.At; r.callsPsi=c.Psi; r.callsPsit=c.Psit;
r.time=op.elapsed(); r.relStep=NaN; r.rankDeficient=0;
end
