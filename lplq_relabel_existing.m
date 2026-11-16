function newFolder = lplq_relabel_existing(runFolder)
% 将一次已有实验复制到新目录，更新输出名称并重新画图，不重新求解。
% 输入runFolder是包含configuration.mat的一次实验根目录，即原outdir。
% 原始目录保留；复制目录内RRE、解、迭代历史、时间等数值不重新计算。
setup_lplq;
runFolder=char(runFolder);
if exist(fullfile(runFolder,'configuration.mat'),'file')~=2
    error('LpLq:RelabelFolder','请选择包含configuration.mat的一次实验目录。');
end
[ok,attr]=fileattrib(runFolder);
if ~ok, error('LpLq:RelabelFolder','无法读取输入目录。'); end
runFolder=attr.Name; [parent,base]=fileparts(runFolder);
[~,suffix]=fileparts(tempname(parent));
newFolder=fullfile(parent,[base,'_MMGKS_labels_',suffix]);
[ok,msg]=copyfile(runFolder,newFolder);
if ~ok, error('LpLq:RelabelCopy','复制失败：%s',msg); end

% CSV只替换完整名称字段，不重新解析/格式化任何数值。
csvFiles=dir(fullfile(newFolder,'**','*.csv'));
for j=1:numel(csvFiles)
    path=fullfile(csvFiles(j).folder,csvFiles(j).name); old=fileread(path);
    changed=strrep(old,'"CAPPED"','"MM-GKS(Kmax)"');
    changed=regexprep(changed,'(?m)(^|,)CAPPED(?=,|\r?$)','$1MM-GKS(Kmax)');
    if ~strcmp(old,changed), write_text(path,changed); end
end

% 更新MAT中的显示字段。cfg.methods和info.options.mode有意保持CAPPED。
matFiles=dir(fullfile(newFolder,'**','results.mat')); plotCount=0;
for j=1:numel(matFiles)
    path=fullfile(matFiles(j).folder,matFiles(j).name); data=load(path);
    if ~isfield(data,'results') || ~iscell(data.results), continue; end
    changed=false;
    for k=1:numel(data.results)
        item=data.results{k};
        if isfield(item,'name')
            label=lplq_method_label(item.name);
            changed=changed || ~strcmp(char(item.name),label); item.name=label;
        end
        if isfield(item,'summary') && isfield(item.summary,'method')
            label=lplq_method_label(item.summary.method);
            changed=changed || ~strcmp(char(item.summary.method),label);
            item.summary.method=label;
        end
        data.results{k}=item;
    end
    if ~changed, continue; end
    save(path,'-struct','data','-v7');
    if all(isfield(data,{'xTrue','data','imageSize'}))
        prob=struct('xTrue',data.xTrue,'d',data.data,'imageSize',data.imageSize);
        % 仅移除复制目录中的旧组合图，避免绘图失败后留下旧标签图。
        names={'reconstructions.png','absolute_errors.png','curves.png','objective.png'};
        for a=1:numel(names)
            imagePath=fullfile(matFiles(j).folder,names{a});
            if exist(imagePath,'file')==2, delete(imagePath); end
        end
        lplq_plot_results(data.results,prob,matFiles(j).folder);
        for a=1:numel(names)
            if exist(fullfile(matFiles(j).folder,names{a}),'file')~=2
                warning('LpLq:RelabelPlot','未生成组合图：%s',names{a});
            end
        end
        plotCount=plotCount+1;
    end
    clear data;
end

% 刷新跨尺寸汇总图；只读取已有数值，不重做重建。
scaleFiles=dir(fullfile(newFolder,'**','scale_summary.csv'));
for j=1:numel(scaleFiles)
    path=fullfile(scaleFiles(j).folder,scaleFiles(j).name);
    T=readtable(path); structs=table2struct(T); rows=cell(numel(structs),1);
    for k=1:numel(structs), rows{k}=structs(k); end
    oldPlot=fullfile(scaleFiles(j).folder,'scaling.png');
    if exist(oldPlot,'file')==2, delete(oldPlot); end
    lplq_plot_scale(rows,scaleFiles(j).folder);
    if exist(oldPlot,'file')~=2
        warning('LpLq:RelabelPlot','未生成跨尺寸汇总图：%s',oldPlot);
    end
end

% 如果原来保存过独立恢复图，只改复制目录里的文件名。
images=dir(fullfile(newFolder,'**','CAPPED_final.png'));
for j=1:numel(images)
    old=fullfile(images(j).folder,images(j).name);
    renamed=fullfile(images(j).folder,'MM-GKS(Kmax)_final.png');
    if exist(renamed,'file')==2
        error('LpLq:RelabelConflict','复制目录中已存在目标文件：%s',renamed);
    end
    [ok,msg]=movefile(old,renamed);
    if ~ok, error('LpLq:RelabelMove','改名失败：%s',msg); end
end
manifests=dir(fullfile(newFolder,'**','run_manifest.txt'));
for j=1:numel(manifests)
    path=fullfile(manifests(j).folder,manifests(j).name);
    write_text(path,strrep(fileread(path),'CAPPED','MM-GKS(Kmax)'));
end
fprintf('显示名称已更新；处理了%d个包含CAPPED的重建结果。\n',plotCount);
fprintf('原始结果保留，新目录：%s\n',newFolder);
end
function write_text(path,content)
% 原文件中的数值字符串保持原样，只写回替换后的文本。
fid=fopen(path,'w');
if fid<0, error('LpLq:RelabelWrite','无法写入：%s',path); end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',content);
end
