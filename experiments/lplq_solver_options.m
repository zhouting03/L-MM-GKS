function o = lplq_solver_options(cfg,lambda,truth)
% 所有同目标方法共用 p,q,lambda,eps,初值和初始化步骤。
o.p=cfg.p; o.q=cfg.q; o.lambda=lambda; o.epsP=cfg.epsP; o.epsQ=cfg.epsQ;
o.kmin=cfg.kmin; o.kmax=cfg.kmax; o.maxIter=cfg.maxIter; o.initDim=cfg.initDim;
o.memoryBudgetMiB=cfg.memoryBudgetMiB; o.qrBlockRows=cfg.qrBlockRows;
o.xTrue=truth; o.x0=zeros(numel(truth),1); o.verbose=false;
end
