# 参数来源与修改位置

求解器默认值在`solver/lplq_defaults.m`，实验profile在`experiments/lplq_experiment_defaults.m`。后者显式覆盖前者的部分值。没有注明来自论文的数值，都是本包为可复现实验设置的起点，不是理论最优值。

| 参数 | 默认与来源 | 含义或修改原则 |
|---|---|---|
| `x0` | 零向量 | 无先验时的共同起点；流式后续块可用上一块解 |
| `V0` | 空 | 自动两阶段初始化；非空时先加入当前问题的x和完整g再正交化 |
| `p` | 1 | 鲁棒数据项；Gaussian可测试2，Laplace可测试1；小于1会非凸 |
| `q` | 1 | 梯度稀疏先验；小于1是非凸备选，不能保证RRE更低 |
| `lambda` | 验证后固定 | 主求解器直接调用的默认值0.03只是起点；实验由独立观测选一次 |
| `lambdaGrid` | 0.001、0.003、0.01、0.03、0.1、0.3 | 验证候选；p不等于1时乘已知噪声RMS的p−1次幂；端点最优需检查网格 |
| `pilotSteps` | quick 25 / paper_small 35 | 选参试算次数，使用LM-TSVD；不是正式迭代中的GCV |
| `epsP,epsQ` | 都为0.001 | 针对[0,1]图像的平滑尺度；与噪声标准差不是同一概念，主迭代中固定 |
| `kmax` | 15 / 25 | 最大空间列数，控制主要缓存；对应原文容量设置思想 |
| `kmin` | 5 | 压缩后至多保留的维数；默认包含x、g和最多3个历史方向 |
| `kminGrid` | quick 3、5、8；paper_small 5、10、15 | 后一组参考论文去模糊测试；均须小于kmax |
| `initDim` | 求解器15；quick实验13 | 初始GKB步数，实际截为min(initDim,kmax−2)；预留x和g |
| `init` | two-stage | 一次GKB代理解后重新建立双权重Krylov种子；simple用于简化对照 |
| `maxIter` | 60 / 150 | 总代理投影求解预算，包含一次初始化代理求解 |
| `gradRelTol,gradAbsTol` | 1e−6、1e−10 | 终止阈值absTol+relTol*max(1,初始完整梯度范数) |
| `orthTol` | 1e−12 | 双重MGS的相关性/新增方向判据，非目标函数容差 |
| `cacheTol` | 1e−10 | 检查新基属于旧空间后才通过线性组合复用算子缓存 |
| `refreshEvery` | 10 | 每10次求解用真实算子刷新当前解的正向像 |
| `descentTol` | 5e−12 | 浮点目标/代理下降检查的相对容差 |
| `qrBlockRows` | 1024 | 分块QR的行数；小块更省临时内存，可能更慢 |
| `memoryBudgetMiB` | 256 | 工作数组估算的拒绝阈值，不是MATLAB进程或操作系统内存上限 |
| `returnBasis` | false | 流式实验才返回至多kmin列的可回收基；默认减少输出体积 |
| `checkAdjoint` | true | 求解前检查A/At、Psi/Psit的内积一致性，调用次数计入统计 |

## 压缩策略参数

| 参数 | 默认 | 解释 |
|---|---|---|
| `compression` | TSVD | 取当前双权重小矩阵最大奇异值对应的右奇异方向 |
| `rbdTol` | 1e−5 | 参考原文RBD容差；本实现按相对最大列残差停止 |
| `rbdSeed` | 17 | 仅控制RBD首列；每轮压缩确定性变化并恢复外部rng |
| `socTol` | 0 | 当前图像归一化后不直接照搬原文阈值1；最大选择kmin−2个分量 |
| `secRho` | 当前固定lambda | 辅助小问题的启发式罚系数，不是额外修改主问题的lambda |
| `secEps` | 1e−4 | SEC辅助L1问题的平滑尺度 |
| `secMaxIter,secTol` | 12、1e−5 | SEC辅助MM预算及相邻系数相对变化容差；只用于选方向 |

只允许已知字段；拼写错误会报错。通过实验入口改配置，例如：

```matlab
c = struct();
c.p = 1; c.q = 0.7;
c.kmin = 5; c.kmax = 20; c.maxIter = 100;
c.lambdaPolicy = 'fixed'; c.lambdaFixed = 0.1;
c.methods = {'FULL','CAPPED','RESTART','TSVD'};
outdir = run_lplq_experiments('compression','quick',c);
```

直接调用求解器时字段名是`lambda`，实验入口的固定选参字段是`lambdaFixed`；二者不要混淆。没有每步选lambda、自动改变epsilon或自动裁剪迭代解的隐藏操作。
