function [x,info] = lplq_lmmgks(op,d,opts)
% 固定参数双权重有限内存 MM-GKS；算子输入只需支持实列向量。
% op.n,m,r,A,At,Psi,Psit；opts 见 lplq_defaults，主循环从不更新 lambda/eps。
if nargin<3, opts=struct(); end
for name={'n','m','r','A','At','Psi','Psit'}
    if ~isfield(op,name{1}), error('LpLq:Operator','Missing op.%s',name{1}); end
end
for name={'n','m','r'}
    validateattributes(op.(name{1}),{'numeric'},{'scalar','integer','positive'});
end
d=double(d(:));
if numel(d)~=op.m || ~isreal(d) || any(~isfinite(d)), error('LpLq:Data','Invalid data vector.'); end
o=lplq_defaults(opts,op.n); plan=lplq_memory_plan(op,o);
clock0=tic; counts=struct('A',0,'At',0,'Psi',0,'Psit',0);
f=op; f.A=@forward; f.At=@adjoint; f.Psi=@regularizer; f.Psit=@regAdjoint;
f.getCounts=@getCounts; f.elapsed=@elapsed;
if o.checkAdjoint
    v=sin((1:op.n)'*sqrt(2)); a=cos((1:op.m)'*sqrt(3)); b=sin((1:op.r)'*sqrt(5));
    av=f.A(v); lv=f.Psi(v);
    errA=abs(a'*av-v'*f.At(a))/max([1,norm(a)*norm(av)]);
    errL=abs(b'*lv-v'*f.Psit(b))/max([1,norm(b)*norm(lv)]);
    if max(errA,errL)>1e-10, error('LpLq:Adjoint','Adjoint check failed.'); end
end
st=lplq_state(f,d,o.x0,o); target=o.gradAbsTol+o.gradRelTol*max(1,st.gnorm);
S=lplq_cache(f,zeros(op.n,0),plan.capacity);
records=cell(o.maxIter+1,1); records{1}=lplq_record(S,st,o,f); N=0;
status='max_iter';
if st.gnorm<=target
    status='converged_initial';
else
    [S,st,N]=lplq_initialize(f,d,o,st,plan.capacity);
    if N==0, records{1}=lplq_record(S,st,o,f);
    else, records{2}=lplq_record(S,st,o,f); end
    while N<o.maxIter && st.gnorm>target
        if S.k>=o.kmax && ~strcmp(o.mode,'FULL')
            if strcmp(o.mode,'CAPPED'), status='memory_limit'; break; end
            S=lplq_compress(S,st,f,d,o);
            records{N+1}=lplq_record(S,st,o,f);
        end
        if strcmp(o.mode,'FULL'), B=o.maxIter-N;
        else, B=min(o.maxIter-N,o.kmax-S.k); end
        [S,st,rr]=lplq_enlarge(S,st,f,d,o,B,target,N);
        for j=1:numel(rr), records{N+j+1}=rr{j}; end
        N=N+numel(rr);
        if o.verbose
            fprintf('%s/%s iter=%d F=%.5e |g|=%.3e dim=%d comp=%d\n', ...
                o.mode,o.compression,N,st.F,st.gnorm,S.k,S.nCompress);
        end
    end
end
% 最终必须用真实算子复核目标与完整梯度，不能只看相邻解变化。
st=lplq_state(f,d,st.x,o);
if st.gnorm<=target, status='converged_gradient'; end
info.basis=[];
if o.returnBasis && S.k>0
    info.basis=lplq_keep_basis([st.x,st.g],S.V(:,1:S.k),o.kmin,o.orthTol);
end
last=lplq_record(S,st,o,f);
if ~isempty(records{N+1}), last.relStep=records{N+1}.relStep; end
records{N+1}=last;
history=struct(); names=fieldnames(last);
for j=1:numel(names)
    a=zeros(N+1,1); for k=1:N+1, a(k)=records{k}.(names{j}); end
    history.(names{j})=a;
end
history.iteration=(0:N)'; history.lambda=o.lambda*ones(N+1,1);
history.epsP=o.epsP*ones(N+1,1); history.epsQ=o.epsQ*ones(N+1,1);
info.history=history; info.status=status; info.converged=st.gnorm<=target;
info.iterations=N; info.objective=st.F; info.gradNorm=st.gnorm;
info.gradTarget=target; info.counts=counts; info.time=toc(clock0);
info.nCompress=S.nCompress; info.nSafeguard=S.nSafeguard;
info.maxPreserveError=S.maxPreserveError; info.memory=plan;
info.peakBasisDim=S.maxDim; info.peakAllocatedCacheMiB=S.peakAllocatedMiB;
o.xTrue=[]; o.x0=[]; o.V0=[]; info.options=o; x=st.x;

    function y=forward(v)
        counts.A=counts.A+1; y=check(op.A(v),op.m);
    end
    function y=adjoint(v)
        counts.At=counts.At+1; y=check(op.At(v),op.n);
    end
    function y=regularizer(v)
        counts.Psi=counts.Psi+1; y=check(op.Psi(v),op.r);
    end
    function y=regAdjoint(v)
        counts.Psit=counts.Psit+1; y=check(op.Psit(v),op.n);
    end
    function c=getCounts(), c=counts; end
    function t=elapsed(), t=toc(clock0); end
end
function y=check(y,n)
y=double(y(:));
if numel(y)~=n || ~isreal(y) || any(~isfinite(y)), error('LpLq:Operator','Invalid operator output.'); end
end
