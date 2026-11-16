function lplq_plot_scale(rows,outdir)
% 跨图像尺寸汇总精度、实际计时、缓存及工作数组估算，避免混称进程峰值。
labels={'FULL',lplq_method_label('CAPPED'),'TSVD'};
fields={'RRE','seconds','peakAllocatedCacheMiB','estimatedWorkMiB'};
ylabels={'Final RRE','Solver elapsed seconds','Peak V/AV/PsiV cache MiB', ...
    'Estimated working arrays MiB'};
try
    fig=figure('Visible','off','Position',[50,50,1050,680]);
    for panel=1:4
        subplot(2,2,panel); hold on;
        for method=1:numel(labels)
            x=[]; y=[];
            for j=1:numel(rows)
                if strcmp(lplq_method_label(rows{j}.method),labels{method})
                    x(end+1)=rows{j}.imageSize^2; y(end+1)=rows{j}.(fields{panel}); %#ok<AGROW>
                end
            end
            [x,idx]=sort(x); plot(x,y(idx),'-o','LineWidth',1.2);
        end
        xlabel('Unknowns (image side squared)'); ylabel(ylabels{panel});
        grid on; legend(labels,'Interpreter','none','Location','best'); hold off;
    end
    print(fig,fullfile(outdir,'scaling.png'),'-dpng','-r140'); close(fig);
catch ME
    warning('LpLq:Plot','Scale figure unavailable: %s',ME.message);
    if exist('fig','var') && ishghandle(fig), close(fig); end
end
end
