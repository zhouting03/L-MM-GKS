function V = lplq_krylov_seed(op,st,d,o,steps)
% 全重正交的对称 Krylov 递推：K(H,b)，不将 H 显式构造成矩阵。
b=op.At(st.wd.*d); nb=norm(b); V=zeros(op.n,steps);
if steps==0 || nb==0, V=zeros(op.n,0); return; end
v=b/nb; used=0;
for j=1:steps
    used=used+1; V(:,used)=v;
    w=op.At(st.wd.*op.A(v))+o.lambda*op.Psit(st.wr.*op.Psi(v));
    Q=V(:,1:used); w=w-Q*(Q'*w); w=w-Q*(Q'*w); nw=norm(w);
    if nw<=o.orthTol, break; end
    v=w/nw;
end
V=V(:,1:used);
end
