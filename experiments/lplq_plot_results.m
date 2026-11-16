function lplq_plot_results(results,prob,outdir)
% 统一色轴重建/误差图；RRE 对次数/时间，缓存与完整梯度曲线。
try
    n=numel(results);
    labels=cellfun(@(r) lplq_method_label(r.name),results,'UniformOutput',false);
    colors=lines(n); styles={'-','--',':','-.','--','-.',':'};
    fig=figure('Visible','off','Position',[50,50,1300,620]);
    images=[{prob.xTrue},{prob.d},cellfun(@(r)r.x,results,'UniformOutput',false)'];
    titles=[{'Truth','Observed'},labels'];
    for j=1:numel(images)
        subplot(2,ceil(numel(images)/2),j); imagesc(reshape(images{j},prob.imageSize),[0,1]);
        axis image off; title(titles{j},'Interpreter','none');
    end
    colormap(gray(256)); print(fig,fullfile(outdir,'reconstructions.png'),'-dpng','-r140'); close(fig);
    maxerr=0; for j=1:n, maxerr=max(maxerr,max(abs(results{j}.x-prob.xTrue))); end
    fig=figure('Visible','off','Position',[50,50,1300,400]);
    for j=1:n
        subplot(1,n,j); imagesc(reshape(abs(results{j}.x-prob.xTrue),prob.imageSize),[0,max(maxerr,eps)]);
        axis image off; title(sprintf('%s RRE=%.3f',labels{j},results{j}.metrics.RRE),'Interpreter','none');
    end
    colormap(gray(256)); print(fig,fullfile(outdir,'absolute_errors.png'),'-dpng','-r140'); close(fig);
    fig=figure('Visible','off','Position',[50,50,1100,720]);
    for k=1:4
        subplot(2,2,k); hold on;
        for j=1:n
            h=results{j}.info.history;
            switch k
                case 1, x=h.iteration; y=max(h.rre,eps);
                case 2, x=h.time; y=max(h.rre,eps);
                case 3, x=h.iteration; y=h.allocatedCacheMiB;
                case 4, x=h.iteration; y=max(h.gradNorm,eps);
            end
            plot(x,y,'Color',colors(j,:),'LineStyle',styles{1+mod(j-1,numel(styles))},'LineWidth',1.3);
            % 终点标记让容量提前停止的短曲线也清晰可见。
            plot(x(end),y(end),'o','Color',colors(j,:),'HandleVisibility','off');
        end
        if k~=3, set(gca,'YScale','log'); end
        grid on; legend(labels,'Interpreter','none','Location','best');
        if k==1, xlabel('Projected solves including initialization'); ylabel('RRE');
        elseif k==2, xlabel('Solver elapsed seconds'); ylabel('RRE');
        elseif k==3, xlabel('Projected solves'); ylabel('Allocated V/AV/PsiV cache MiB');
        else, xlabel('Projected solves'); ylabel('Complete gradient norm'); end
        hold off;
    end
    print(fig,fullfile(outdir,'curves.png'),'-dpng','-r140'); close(fig);
    fig=figure('Visible','off'); hold on;
    for j=1:n
        h=results{j}.info.history;
        plot(h.iteration,h.objective,'Color',colors(j,:), ...
            'LineStyle',styles{1+mod(j-1,numel(styles))},'LineWidth',1.3);
    end
    set(gca,'YScale','log');
    xlabel('Projected solves'); ylabel('Same fixed objective'); grid on;
    legend(labels,'Interpreter','none'); print(fig,fullfile(outdir,'objective.png'),'-dpng','-r140'); close(fig);
catch ME
    warning('LpLq:Plot','Figures unavailable, numeric CSV/MAT already saved: %s',ME.message);
    if exist('fig','var') && ishghandle(fig), close(fig); end
end
end
