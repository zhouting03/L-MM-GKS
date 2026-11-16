function V = lplq_gkb_seed(op,d,wd,steps,tol)
% 对 D*A、D*d 做 GKB，返回右基；种子为零或发生 breakdown 时正常截断。
V=zeros(op.n,steps); U=zeros(op.m,steps+1);
sw=sqrt(wd); f=sw.*d; beta=norm(f);
if steps==0 || beta==0, V=zeros(op.n,0); return; end
u=f/beta; U(:,1)=u; oldv=zeros(op.n,1); used=0;
for j=1:steps
    a=op.At(sw.*u)-beta*oldv;
    if used>0
        Q=V(:,1:used); a=a-Q*(Q'*a); a=a-Q*(Q'*a);
    end
    alpha=norm(a);
    if alpha<=tol, break; end
    v=a/alpha; used=used+1; V(:,used)=v;
    b=sw.*op.A(v)-alpha*u;
    Q=U(:,1:j); b=b-Q*(Q'*b); b=b-Q*(Q'*b);
    beta=norm(b);
    if beta<=tol, break; end
    oldv=v; u=b/beta; U(:,j+1)=u;
end
V=V(:,1:used);
end
