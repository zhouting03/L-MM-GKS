function plan = lplq_memory_plan(op,o)
% 工作数组估算：不是 MATLAB 进程峰值，也不是操作系统强制内存上限。
if strcmp(o.mode,'FULL')
    cap=min(op.n,max(o.kmax,o.kmin+o.maxIter));
else
    cap=o.kmax;
end
n=double(op.n); m=double(op.m); r=double(op.r);
b=min(double(o.qrBlockRows),max(m,r));
plan.capacity=cap;
plan.cacheAtCapacityMiB=8*(n+m+r)*cap/2^20;
% 预留增长/压缩副本、分块 QR 工作区及若干完整向量。
plan.estimatedWorkMiB=8*(2*(n+m+r)*cap+6*b*(cap+1)+ ...
    24*(n+m+r)+24*(cap+1)^2)/2^20;
plan.budgetMiB=o.memoryBudgetMiB;
if plan.estimatedWorkMiB>o.memoryBudgetMiB
    error('LpLq:MemoryBudget', ...
        'Estimated working arrays %.1f MiB exceed budget %.1f MiB. Reduce size/iterations/kmax.', ...
        plan.estimatedWorkMiB,o.memoryBudgetMiB);
end
end
