function lplq_streaming_demo(cfg,outdir)
% 参考算法4.1：顺序处理当前数据块，回收解和空间；不宣称等价于全数据最小化。
n=min(cfg.imageSize,48); truth=lplq_phantom(n); truth=truth(:);
angles={0:6:42,48:6:90,96:12:180}; ng=numel(angles);
Xwarm=[]; Vwarm=[]; rows={}; results=cell(3,ng);
for group=1:ng
    [A,~]=lplq_tomo_matrix(n,angles{group});
    op=make_op(A,n); clean=op.A(truth);
    [data,noise]=lplq_noise(clean,'gaussian',0.005,cfg.seed+group,0);
    [validation,~]=lplq_noise(clean,'gaussian',0.005,cfg.seed+group+1009,0);
    c=cfg; c.p=2; c.q=1; c.noiseType='gaussian'; c.pilotSteps=min(25,cfg.pilotSteps);
    c.maxIter=min(cfg.maxIter,80); c.kmin=min(cfg.kmin,5); c.kmax=min(cfg.kmax,15);
    prob=struct('op',op,'d',data,'validationData',validation,'noise',noise,'xTrue',truth);
    [lambda,cal]=lplq_select_lambda(prob,c); %#ok<NASGU>
    for method=1:3
        o=lplq_solver_options(c,lambda,truth); o.returnBasis=true; o.V0=[];
        if method==1
            label='independent'; o.x0=zeros(numel(truth),1); o.V0=[];
        elseif method==2
            label='warm_solution';
            if group>1, o.x0=Xwarm{2}; end
            o.V0=[];
        else
            label='recycle_solution_basis';
            if group>1, o.x0=Xwarm{3}; o.V0=Vwarm; end
        end
        [x,info]=lplq_lmmgks(op,data,o);
        if isempty(Xwarm), Xwarm=cell(3,1); end
        Xwarm{method}=x; if method==3, Vwarm=info.basis; end
        row.method=label; row.group=group; row.RRE=norm(x-truth)/norm(truth);
        row.lambda=lambda; row.seconds=info.time; row.solves=info.iterations;
        row.cacheMiB=info.peakAllocatedCacheMiB; row.recycledDim=size(o.V0,2);
        rows{end+1}=row; %#ok<AGROW>
        results{method,group}=struct('x',x,'info',info,'summary',row);
        fprintf('CT group %d %-24s RRE %.4f\n',group,label,row.RRE);
    end
    clear A op prob; % 本块算子不传入结果，下一块重新构造。
end
% 流式阶段结束之后，才重新组装小规模全数据参考。
allA=cell(ng,1); alld=cell(ng,1); allval=cell(ng,1); noiseEnergy=0; count=0;
for group=1:ng
    allA{group}=lplq_tomo_matrix(n,angles{group}); clean=allA{group}*truth;
    [alld{group},ni]=lplq_noise(clean,'gaussian',0.005,cfg.seed+group,0);
    allval{group}=lplq_noise(clean,'gaussian',0.005,cfg.seed+group+1009,0);
    noiseEnergy=noiseEnergy+ni.rms^2*numel(clean); count=count+numel(clean);
end
A=vertcat(allA{:}); data=vertcat(alld{:}); validation=vertcat(allval{:}); op=make_op(A,n);
prob=struct('op',op,'d',data,'validationData',validation,'noise',struct('rms',sqrt(noiseEnergy/count)));
[lambdaBatch,batchCalibration]=lplq_select_lambda(prob,c); %#ok<NASGU>
o=lplq_solver_options(c,lambdaBatch,truth); o.returnBasis=false;
[xBatch,infoBatch]=lplq_lmmgks(op,data,o); %#ok<NASGU>
row.method='batch_reference'; row.group=0; row.RRE=norm(xBatch-truth)/norm(truth);
row.lambda=lambdaBatch; row.seconds=infoBatch.time; row.solves=infoBatch.iterations;
row.cacheMiB=infoBatch.peakAllocatedCacheMiB; row.recycledDim=0; rows{end+1}=row;
lplq_write_csv(fullfile(outdir,'streaming_summary.csv'),rows);
save(fullfile(outdir,'streaming.mat'),'results','truth','angles','xBatch','infoBatch','batchCalibration','-v7');
if cfg.makeFigures
    try
        fig=figure('Visible','off','Position',[50,50,1050,650]);
        shortLabels={'Independent','Warm x','Recycle x,V'};
        for g=1:ng
            for j=1:3
                subplot(3,ng,(j-1)*ng+g); imagesc(reshape(results{j,g}.x,n,n),[0,1]); axis image off;
                title(sprintf('%s | block%d | RRE %.3f',shortLabels{j},g,results{j,g}.summary.RRE),'Interpreter','none');
            end
        end
        colormap(gray(256)); print(fig,fullfile(outdir,'streaming_reconstructions.png'),'-dpng','-r140'); close(fig);
    catch ME
        warning('LpLq:Plot','Streaming plot unavailable: %s',ME.message);
        if exist('fig','var') && ishghandle(fig), close(fig); end
    end
end
end
function op=make_op(A,n)
op.n=n*n; op.m=size(A,1); op.r=2*n*n;
op.A=@(x) A*x; op.At=@(y) A'*y;
op.Psi=@(x) lplq_grad2(x,n,n,false); op.Psit=@(y) lplq_grad2(y,n,n,true);
end
