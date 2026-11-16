function [x,info] = lplq_mmgks(op,d,opts)
% 同目标、同初始化、同新梯度的无压缩基线；不是原文 GCV 软件的逐字复现。
if nargin<3, opts=struct(); end
opts.mode='FULL'; [x,info]=lplq_lmmgks(op,d,opts);
end
