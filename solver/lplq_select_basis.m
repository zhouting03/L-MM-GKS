function W = lplq_select_basis(fac,a,t,o,cycle)
% 四种真实压缩策略：RBD 缺失时不冒充 TSVD；SEC 仅用于选方向。
k=size(fac.H,2);
if t==0, W=zeros(k,0); return; end
t=min(t,k);
switch o.compression
    case 'TSVD'
        [~,~,U]=svd(fac.H,'econ'); W=U(:,1:t);
    case 'RBD'
        X=fac.H'; W=zeros(k,0); R=X;
        base=max(sqrt(sum(X.^2,1))); pool=find(sum(X.^2,1)>0);
        if isempty(pool), return; end
        state=rng; rng(o.rbdSeed+cycle,'twister'); j=pool(randi(numel(pool))); rng(state);
        for it=1:t
            v=R(:,j);
            if ~isempty(W), v=v-W*(W'*v); v=v-W*(W'*v); end
            nv=norm(v); if nv<=o.orthTol*max(base,realmin), break; end
            W=[W,v/nv]; %#ok<AGROW>
            R=X-W*(W'*X); norms=sqrt(sum(R.^2,1)); [val,j]=max(norms);
            if val<=o.rbdTol*max(base,realmin), break; end
        end
    case 'SOC'
        W=coordinate_selection(a,t,o.socTol,k);
    case 'SEC'
        z=a;
        for it=1:o.secMaxIter
            rz=fac.Rr*z; weights=hypot(rz,o.secEps).^(-1);
            H=[fac.Rd;sqrt(o.secRho)*bsxfun(@times,fac.Rr,sqrt(weights))];
            rhs=[fac.cd;zeros(size(fac.Rr,1),1)];
            znew=H\rhs;
            if norm(znew-z)<=o.secTol*max(1,norm(z)), z=znew; break; end
            z=znew;
        end
        W=coordinate_selection(z,t,o.socTol,k);
end
end
function W=coordinate_selection(a,t,tol,k)
[~,idx]=sort(abs(a),'descend'); idx=idx(1:min(t,numel(idx)));
idx=idx(abs(a(idx))>tol); E=eye(k); W=E(:,idx);
end
