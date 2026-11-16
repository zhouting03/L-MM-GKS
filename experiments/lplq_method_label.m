function label = lplq_method_label(method)
% 统一输出名称；求解器和cfg.methods继续使用内部标识CAPPED。
if iscell(method) && isscalar(method), method=method{1}; end
label=char(method);
if strcmpi(strtrim(label),'CAPPED')
    label='MM-GKS(Kmax)';
end
end
