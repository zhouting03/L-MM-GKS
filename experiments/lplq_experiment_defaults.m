function cfg = lplq_experiment_defaults(profile)
% quick 是默认小规模入口；paper_small 增大预算，仍不复现原文大尺寸。
if nargin<1, profile='quick'; end
cfg.profile=char(profile); cfg.imageSize=64; cfg.imageSource='phantom';
cfg.hstPath='/Users/zhouting/Documents/编程/文献代码库/L-MM-GKS/LM-MM-GKS/utilities/HSTgray.jpg';
cfg.allowLargerImages=false; cfg.noiseType='laplace_impulse'; cfg.noiseLevel=0.02;
cfg.impulseFraction=0.02; cfg.seed=17; cfg.p=1; cfg.q=1;
cfg.epsP=1e-3; cfg.epsQ=1e-3; cfg.lambdaPolicy='validation'; cfg.lambdaFixed=0.03;
cfg.lambdaGrid=[1e-3,3e-3,1e-2,3e-2,1e-1,3e-1]; cfg.pilotSteps=25;
%cfg.maxIter=60; 
cfg.maxIter=100; 
cfg.kmin=5; cfg.kmax=15; cfg.initDim=13;
cfg.memoryBudgetMiB=256; cfg.qrBlockRows=1024;
cfg.makeFigures=true; cfg.saveFinalImages=true;
cfg.methods={'FULL','CAPPED','RESTART','TSVD','RBD','SOC','SEC'};
cfg.kminGrid=[3,5,8]; cfg.sizeGrid=[64,96,128]; cfg.repetitions=1;
cfg.modelPairs=[2,1;1,1;1,0.7]; cfg.noiseTypes={'gaussian','laplace','laplace_impulse'};
if strcmp(profile,'paper_small')
    cfg.imageSize=96; cfg.maxIter=150; cfg.kmax=25; cfg.initDim=15;
    cfg.kminGrid=[5,10,15]; cfg.pilotSteps=35;
elseif ~strcmp(profile,'quick')
    error('LpLq:Profile','Use quick or paper_small.');
end
end
