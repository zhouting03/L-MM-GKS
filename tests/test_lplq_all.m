function report = test_lplq_all()
% 实质数值检查：导数、QR、压缩保持、凸参考解、退化情况与内存预检。
setup_lplq(); saved=rng; cleanup=onCleanup(@()rng(saved)); %#ok<NASGU>
rng(23,'twister'); A=randn(31,20)/sqrt(31); L=diff(eye(20));
op=lplq_from_matrices(A,L); d=randn(31,1); x=randn(20,1);
opts=struct('p',0.8,'q',0.7,'lambda',0.07,'epsP',0.03,'epsQ',0.04, ...
    'kmin',4,'kmax',9,'maxIter',24,'checkAdjoint',true,'returnBasis',true);
o=lplq_defaults(opts,op.n); st=lplq_state(op,d,x,o);
v=randn(20,1); v=v/norm(v); h=1e-6;
f1=lplq_state(op,d,x+h*v,o); f2=lplq_state(op,d,x-h*v,o);
assert(abs((f1.F-f2.F)/(2*h)-v'*st.g)<1e-5,'Gradient finite difference failed.');
V=lplq_keep_basis([x,st.g],randn(20,5),7,1e-12); S=lplq_cache(op,V,10);
[z,~]=lplq_projected_qr(S,st,d,o);
H=[bsxfun(@times,A*V,sqrt(st.wd));sqrt(o.lambda)*bsxfun(@times,L*V,sqrt(st.wr))];
rhs=[sqrt(st.wd).*d;zeros(size(L,1),1)]; reference=H\rhs;
assert(norm(z-reference)<1e-8*max(1,norm(reference)),'Block QR disagrees with direct weighted LS.');
pairs=[2,2;2,1;1,1;1,0.7;0.8,0.6]; methods={'FULL','CAPPED','RESTART','TSVD','RBD','SOC','SEC'};
solves=0; maxIncrease=0; maxPreserve=0;
for p=1:size(pairs,1)
    for j=1:numel(methods)
        oo=opts; oo.p=pairs(p,1); oo.q=pairs(p,2);
        if j<=3, oo.mode=methods{j}; else, oo.mode='LM'; oo.compression=methods{j}; end
        [xx,info]=lplq_lmmgks(op,d,oo); solves=solves+1;
        F=info.history.objective; inc=max([0;diff(F)]); maxIncrease=max(maxIncrease,inc);
        assert(inc<1e-8*max(1,F(1)),'Fixed objective increased.');
        assert(all(info.history.lambda==oo.lambda),'lambda changed.');
        assert(all(info.history.epsP==oo.epsP),'epsP changed.');
        assert(info.maxPreserveError<1e-8,'Compression did not preserve x,g.');
        maxPreserve=max(maxPreserve,info.maxPreserveError);
        if ~strcmp(oo.mode,'FULL'), assert(info.peakBasisDim<=oo.kmax,'Capacity exceeded.'); end
        assert(all(isfinite(xx)),'Nonfinite output.');
    end
end
% 充分迭代的二次模型与精确 Tikhonov 解比较。
oo=opts; oo.p=2; oo.q=2; oo.mode='FULL'; oo.maxIter=40;
oo.gradRelTol=1e-11; oo.gradAbsTol=1e-12;
[xx,info]=lplq_lmmgks(op,d,oo); exact=(A'*A+oo.lambda*L'*L)\(A'*d);
quadraticError=norm(xx-exact)/max(1,norm(exact));
assert(quadraticError<1e-7,'Quadratic reference solution mismatch.');
% kmin=2 时 tSVD 无历史方向，应与两方向 restart 一致。
oo=opts; oo.p=1; oo.q=1; oo.kmin=2; oo.mode='LM'; oo.compression='TSVD';
[x1,i1]=lplq_lmmgks(op,d,oo); oo.mode='RESTART'; [x2,i2]=lplq_lmmgks(op,d,oo);
assert(norm(x1-x2)<1e-8 && i1.iterations==i2.iterations,'Two-direction equivalence failed.');
% 零数据零初值应直接驻定，不能 GKB 除零。
[xx,info]=lplq_lmmgks(op,zeros(op.m,1),opts);
assert(norm(xx)==0 && info.iterations==0 && info.converged,'Zero-data edge case failed.');
% 长宽不等且每块仅1行，检查矩形 QR。
oo=opts; oo.qrBlockRows=1; oo.maxIter=4; [~,info]=lplq_lmmgks(op,d,oo);
assert(all(isfinite(info.history.objective)),'One-row QR failed.');
% 内存限制必须在分配大型缓存前拒绝。
oo=opts; oo.memoryBudgetMiB=1e-6; caught=false;
try, lplq_lmmgks(op,d,oo); catch ME, caught=strcmp(ME.identifier,'LpLq:MemoryBudget'); end
assert(caught,'Memory guard did not reject the request.');
% 伴随算子对 FFT 与周期差分的检查由实际求解器执行。
cfg=lplq_experiment_defaults('quick'); cfg.imageSize=24;
prob=lplq_make_deblur(cfg); oo=lplq_solver_options(cfg,0.03,prob.xTrue); oo.maxIter=8;
[~,info]=lplq_lmmgks(prob.op,prob.d,oo);
assert(all(diff(info.history.objective)<1e-7),'FFT deblurring descent failed.');
report=struct('passed',true,'modelMethodRuns',solves,'maxObjectiveIncrease',maxIncrease, ...
    'maxPreserveError',maxPreserve,'quadraticRelativeError',quadraticError, ...
    'runtime',version);
fprintf('PASS: %d model/method runs; max objective increase %.3e; preserve error %.3e\n', ...
    solves,maxIncrease,maxPreserve);
end
