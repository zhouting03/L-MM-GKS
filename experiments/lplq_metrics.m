function met = lplq_metrics(x,truth,sz)
% RRE/PSNR 使用原始解；显示图的 [0,1] 色轴不改变求解结果。
x=x(:); truth=truth(:); mse=mean((x-truth).^2);
met.RRE=norm(x-truth)/max(norm(truth),realmin);
met.PSNR=10*log10(1/max(mse,realmin)); met.SSIM=NaN; met.HaarPSI=NaN;
met.SSIMSource='unavailable'; met.HaarPSISource='unavailable';
X=reshape(x,sz); T=reshape(truth,sz);
if exist('ssim','file')==2
    try
        met.SSIM=ssim(X,T,'DynamicRange',1); met.SSIMSource='MATLAB ssim';
    catch
        % 没有兼容的 SSIM 实现时保留 NaN，不以其他分数冒名替代。
    end
end
if exist('HaarPSI','file')==2
    try
        hp=HaarPSI(255*T,255*min(max(X,0),1));
        if isscalar(hp) && isfinite(hp)
            met.HaarPSI=hp; met.HaarPSISource='external HaarPSI; 2D; display clipped to [0,255]';
        end
    catch
        % 外部接口不匹配时不输出可能失真的接近 1 的分数。
    end
end
end
