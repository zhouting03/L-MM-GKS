function [next,meta] = lplq_step(S,st,op,d,o)
% 一次当前代理投影求解；检测数值失败时比较真实梯度 Cauchy 点。
[z,fac]=lplq_projected_qr(S,st,d,o);
ax=S.Y(:,1:S.k)*z; lx=S.Z(:,1:S.k)*z; x=S.V(:,1:S.k)*z;
q0=0.5*sum(st.wd.*st.e.^2)+0.5*o.lambda*sum(st.wr.*st.u.^2);
q1=0.5*sum(st.wd.*(ax-d).^2)+0.5*o.lambda*sum(st.wr.*lx.^2);
next=lplq_state(op,d,x,o,ax,lx); meta.safeguard=false;
if fac.rankDeficient || q1>q0+o.descentTol*max(1,abs(q0)) || ...
        next.F>st.F+o.descentTol*max(1,abs(st.F))
    ag=op.A(st.g); lg=op.Psi(st.g);
    denom=sum(st.wd.*ag.^2)+o.lambda*sum(st.wr.*lg.^2);
    if denom>0 && st.gnorm>0
        alpha=st.gnorm^2/denom;
        cx=st.x-alpha*st.g; ca=st.ax-alpha*ag; cl=st.u-alpha*lg;
        qc=0.5*sum(st.wd.*(ca-d).^2)+0.5*o.lambda*sum(st.wr.*cl.^2);
        if qc<q1 || ~isfinite(q1)
            next=lplq_state(op,d,cx,o,ca,cl); meta.safeguard=true;
        end
    end
end
if next.F>st.F+10*o.descentTol*max(1,abs(st.F))
    error('LpLq:Descent','Surrogate/gradient consistency failure: objective increased.');
end
meta.relStep=norm(next.x-st.x)/max(1,norm(st.x));
meta.rankDeficient=fac.rankDeficient;
end
