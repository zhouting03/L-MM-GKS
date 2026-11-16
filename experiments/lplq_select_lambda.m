function [lambda,detail] = lplq_select_lambda(prob,cfg)
% 两次独立观测的验证选参；不使用真实图像计算选参分数。选完即固定。
detail.grid=cfg.lambdaGrid; detail.score=NaN(size(detail.grid));
detail.totalPilotSeconds=0; detail.policy=cfg.lambdaPolicy;
if strcmp(cfg.lambdaPolicy,'fixed')
    lambda=cfg.lambdaFixed; detail.selected=lambda;
    detail.grid=lambda; detail.score=NaN; detail.selectedIndex=1; detail.onBoundary=false; return;
end
if ~strcmp(cfg.lambdaPolicy,'validation') || isempty(prob.validationData)
    error('LpLq:Calibration','Validation policy requires an independent observation vector.');
end
% p=2 的 lambda 量级通常比 p=1 小；这只调整候选网格，不更新正式迭代参数。
scale=max(prob.noise.rms,1e-3)^(cfg.p-1); detail.grid=cfg.lambdaGrid*scale;
for j=1:numel(detail.grid)
    o=lplq_solver_options(cfg,detail.grid(j),[]); o.x0=zeros(prob.op.n,1);
    o.xTrue=[]; o.maxIter=cfg.pilotSteps; o.mode='LM'; o.compression='TSVD';
    [x,info]=lplq_lmmgks(prob.op,prob.d,o);
    residual=prob.op.A(x)-prob.validationData;
    if strcmp(cfg.noiseType,'gaussian'), detail.score(j)=mean(residual.^2);
    else, detail.score(j)=mean(abs(residual)); end
    detail.totalPilotSeconds=detail.totalPilotSeconds+info.time;
end
[~,idx]=min(detail.score); lambda=detail.grid(idx);
detail.selected=lambda; detail.selectedIndex=idx;
detail.onBoundary=idx==1 || idx==numel(detail.grid);
if detail.onBoundary
    fprintf('  Calibration optimum lies on candidate-grid boundary; inspect lambda_grid.csv.\n');
end
end
