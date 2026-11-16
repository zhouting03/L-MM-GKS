function results = lplq_compare(prob,cfg,lambda,outdir)
% 同目标比较：FULL 与 LM 用同一总预算；CAPPED 单列为存储上限对照。
if ~exist(outdir,'dir'), mkdir(outdir); end
results=cell(numel(cfg.methods),1); rows=results;
for j=1:numel(cfg.methods)
    name=cfg.methods{j}; displayName=lplq_method_label(name);
    o=lplq_solver_options(cfg,lambda,prob.xTrue);
    if any(strcmp(name,{'FULL','CAPPED','RESTART'})), o.mode=name;
    else, o.mode='LM'; o.compression=name; end
    [x,info]=lplq_lmmgks(prob.op,prob.d,o);
    met=lplq_metrics(x,prob.xTrue,prob.imageSize);
    if cfg.saveFinalImages
        imwrite(min(max(reshape(x,prob.imageSize),0),1),fullfile(outdir,[displayName,'_final.png']));
    end
    row.method=displayName; row.p=cfg.p; row.q=cfg.q; row.lambda=lambda;
    row.epsP=cfg.epsP; row.epsQ=cfg.epsQ; row.RRE=met.RRE;
    row.PSNR=met.PSNR; row.SSIM=met.SSIM; row.HaarPSI=met.HaarPSI;
    row.iterations=info.iterations; row.status=info.status; row.seconds=info.time;
    row.objective=info.objective; row.gradNorm=info.gradNorm;
    row.peakBasisDim=info.peakBasisDim; row.nCompress=info.nCompress;
    row.peakAllocatedCacheMiB=info.peakAllocatedCacheMiB;
    row.estimatedWorkMiB=info.memory.estimatedWorkMiB;
    row.callsA=info.counts.A; row.callsAt=info.counts.At;
    row.callsPsi=info.counts.Psi; row.callsPsit=info.counts.Psit;
    row.preserveError=info.maxPreserveError; row.nSafeguard=info.nSafeguard;
    row.SSIMSource=met.SSIMSource; row.HaarPSISource=met.HaarPSISource;
    rows{j}=row; results{j}=struct('name',displayName,'x',x,'info',info,'metrics',met,'summary',row);
    fprintf('%-16s RRE=%.4f PSNR=%5.2f iter=%3d dim=%3d cache=%6.2fMiB time=%6.2fs %s\n', ...
        displayName,met.RRE,met.PSNR,info.iterations,info.peakBasisDim, ...
        info.peakAllocatedCacheMiB,info.time,info.status);
end
lplq_write_csv(fullfile(outdir,'summary.csv'),rows);
% 保存数据与最终解，不保存带函数句柄的 op 或每步完整图像。
xTrue=prob.xTrue; data=prob.d; imageSize=prob.imageSize; noise=prob.noise; source=prob.imageSource; %#ok<NASGU>
save(fullfile(outdir,'results.mat'),'results','cfg','lambda','xTrue','data','imageSize','noise','source','-v7');
if cfg.makeFigures, lplq_plot_results(results,prob,outdir); end
end
