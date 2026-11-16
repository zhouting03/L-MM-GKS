function outdir = run_lplq_experiments(suite,profile,overrides)
% 入口示例：run_lplq_experiments('compression','quick')。
% suite: compression / kmin / memory / robustness / streaming / all。
root=setup_lplq();
if nargin<1 || isempty(suite), suite='compression'; end
if nargin<2 || isempty(profile), profile='quick'; end
if nargin<3, overrides=struct(); end
cfg=lplq_experiment_defaults(profile); names=fieldnames(overrides);
for j=1:numel(names)
    if ~isfield(cfg,names{j}), error('LpLq:Config','Unknown configuration: %s',names{j}); end
    cfg.(names{j})=overrides.(names{j});
end
validateattributes(cfg.imageSize,{'numeric'},{'scalar','integer','>=',16});
if ~any(strcmp(suite,{'compression','kmin','memory','robustness','streaming','all'}))
    error('LpLq:Suite','Unknown experiment suite.');
end
base=fullfile(root,'results'); if ~exist(base,'dir'), mkdir(base); end
% tempname 提供唯一后缀；重复执行不覆盖已有结果。
[~,suffix]=fileparts(tempname(base));
outdir=fullfile(base,[datestr(now,'yyyymmdd_HHMMSS'),'_',char(suite),'_',suffix]); mkdir(outdir);
fid=fopen(fullfile(outdir,'run_manifest.txt'),'w');
fprintf(fid,'Suite: %s\nProfile: %s\nMATLAB/Octave: %s\n',suite,profile,version);
fprintf(fid,'All main iterations use fixed p,q,lambda,epsP,epsQ.\n');
fprintf(fid,'FULL and LM: same solve budget. MM-GKS(Kmax): same basis cap, fewer solves.\n');
fprintf(fid,'Independent-observation validation calibration is separate from solver timing.\n'); fclose(fid);
save(fullfile(outdir,'configuration.mat'),'cfg','-v7');
if strcmp(suite,'all'), suites={'compression','kmin','memory','robustness','streaming'};
else, suites={suite}; end
for s=1:numel(suites)
    name=suites{s}; sub=fullfile(outdir,name); mkdir(sub);
    fprintf('\n=== %s / %s ===\n',name,profile);
    switch name
        case 'compression'
            prob=lplq_make_deblur(cfg); fprintf('Image: %s\n',prob.imageSource);
            [lambda,cal]=lplq_select_lambda(prob,cfg); save_calibration(cal,sub);
            fprintf('Fixed lambda: %.3g; calibration %.2fs (not charged to solver rows).\n',lambda,cal.totalPilotSeconds);
            lplq_compare(prob,cfg,lambda,sub);
        case 'kmin'
            prob=lplq_make_deblur(cfg); [lambda,cal]=lplq_select_lambda(prob,cfg); save_calibration(cal,sub);
            for j=1:numel(cfg.kminGrid)
                c=cfg; c.kmin=cfg.kminGrid(j);
                if c.kmin>=c.kmax, error('LpLq:Config','Each kmin must be smaller than kmax.'); end
                lplq_compare(prob,c,lambda,fullfile(sub,sprintf('kmin_%d',c.kmin)));
            end
        case 'memory'
            rows={};
            for j=1:numel(cfg.sizeGrid)
                c=cfg; c.imageSize=cfg.sizeGrid(j); c.methods={'FULL','CAPPED','TSVD'};
                prob=lplq_make_deblur(c); [lambda,cal]=lplq_select_lambda(prob,c);
                caseDir=fullfile(sub,sprintf('n_%d',c.imageSize)); mkdir(caseDir); save_calibration(cal,caseDir);
                rr=lplq_compare(prob,c,lambda,caseDir);
                for a=1:numel(rr), row=rr{a}.summary; row.imageSize=c.imageSize; rows{end+1}=row; end %#ok<AGROW>
                clear rr prob;
            end
            lplq_write_csv(fullfile(sub,'scale_summary.csv'),rows);
            if cfg.makeFigures, lplq_plot_scale(rows,sub); end
        case 'robustness'
            rows={};
            for j=1:numel(cfg.noiseTypes)
                for rep=1:cfg.repetitions
                    c=cfg; c.noiseType=cfg.noiseTypes{j}; c.seed=cfg.seed+rep-1;
                    prob=lplq_make_deblur(c); % 同一噪声样本用于所有 p,q 模型。
                    for a=1:size(cfg.modelPairs,1)
                        c.p=cfg.modelPairs(a,1); c.q=cfg.modelPairs(a,2); c.methods={'FULL','TSVD'};
                        caseDir=fullfile(sub,sprintf('%s_rep%d_p%g_q%g',c.noiseType,rep,c.p,c.q)); mkdir(caseDir);
                        [lambda,cal]=lplq_select_lambda(prob,c); save_calibration(cal,caseDir);
                        rr=lplq_compare(prob,c,lambda,caseDir);
                        for b=1:numel(rr)
                            row=rr{b}.summary; row.noise=c.noiseType; row.rep=rep; rows{end+1}=row; %#ok<AGROW>
                        end
                    end
                end
            end
            % 不将不同模型的目标值画在同一收敛图中。
            lplq_write_csv(fullfile(sub,'model_summary.csv'),rows);
        case 'streaming'
            lplq_streaming_demo(cfg,sub);
    end
end
fprintf('\nResults: %s\n',outdir);
end
function save_calibration(cal,outdir)
rows=cell(numel(cal.grid),1);
for j=1:numel(cal.grid)
    rows{j}=struct('lambda',cal.grid(j),'validationScore',cal.score(j));
end
lplq_write_csv(fullfile(outdir,'lambda_grid.csv'),rows);
save(fullfile(outdir,'lambda_calibration.mat'),'cal','-v7');
end
