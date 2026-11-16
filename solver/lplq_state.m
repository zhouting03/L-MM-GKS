function st = lplq_state(op,d,x,o,ax,lx)
% 用最新两组权重计算目标与完整梯度；hypot 避免不必要的平方溢出。
if nargin<5, ax=op.A(x); end
if nargin<6, lx=op.Psi(x); end
st.x=x(:); st.ax=ax(:); st.u=lx(:); st.e=st.ax-d;
hp=hypot(st.e,o.epsP); hq=hypot(st.u,o.epsQ);
st.wd=hp.^(o.p-2); st.wr=hq.^(o.q-2);
st.F=sum(hp.^o.p)/o.p+o.lambda*sum(hq.^o.q)/o.q;
st.g=op.At(st.wd.*st.e)+o.lambda*op.Psit(st.wr.*st.u);
st.g=st.g(:); st.gnorm=norm(st.g);
if ~isfinite(st.F) || any(~isfinite(st.g)) || any(~isfinite(st.wd)) || any(~isfinite(st.wr))
    error('LpLq:NonFinite','Nonfinite model state: check scaling and positive smoothing parameters.');
end
end
