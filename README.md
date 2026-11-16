# 双动态权重 Lp–Lq 有限内存 MM–GKS MATLAB 实验包

本包实现前面确认的**固定参数**流程，并提供可直接运行的低内存实验。算法求解

\[
F(x)=\frac1p\sum_i((Ax-d)_i^2+\epsilon_p^2)^{p/2}
+\frac\lambda q\sum_j((\Psi x)_j^2+\epsilon_q^2)^{q/2},\quad 0<p,q\le2.
\]

`lambda, epsP, epsQ, p, q` 在一次正式求解中全部固定。每步更新数据权重和正则权重，使用**最新两组权重下的完整梯度**扩展。压缩直接优先保留高维 `x,g`，与原文“保留解和完整梯度”的思路一致；不是只保留可能遗漏空间外分量的 `V'*g`。

## 1 最快运行

解压后将 MATLAB 当前目录切换到本包根目录：

```matlab
setup_lplq;
report = test_lplq_all();
outdir = run_lplq_experiments('compression','quick');
```

默认 `64×64`、内置 phantom、矩阵自由运动模糊、60次总投影求解、`kmin=5,kmax=15`，按顺序运行7种方法。无需 AIRToolsII、IRTools、HyBRrecycle 或 Optimization Toolbox。默认不需要外部图像。MATLAB 建议 R2020b 或以上；数值核心使用基础 MATLAB 函数。

下一步可运行：

```matlab
run_lplq_experiments('kmin','quick');
run_lplq_experiments('memory','quick');
run_lplq_experiments('robustness','quick');
run_lplq_experiments('streaming','quick');
```

`'all'` 会按顺序执行压缩、保留维数、规模、鲁棒性、流式CT，运行时间明显长于第一条命令，不会同时并行占用多份大型工作区。

## 2 各实验回答什么问题

| 入口 | 对比内容 | 主要用途 |
|---|---|---|
| `compression` | FULL、CAPPED、RESTART、TSVD、RBD、SOC、SEC | 同目标的压缩对比；同时展示限存储基线 |
| `kmin` | `kminGrid` 中不同保留维数 | 对应原文去模糊的保留维数试验 |
| `memory` | 默认64、96、128；FULL/CAPPED/TSVD | RRE、算子次数、缓存及工作数组估算随规模变化 |
| `robustness` | Gaussian/Laplace/Laplace+impulse；(2,1)、(1,1)、(1,0.7) | 区分数据模型效果和优化算法效果 |
| `streaming` | 每块独立、仅热启动解、回收解与空间，另有全数据参考 | 参考原文算法4.1的子空间迁移，默认不超过48×48 |

原文还包含大规模动态PAT。本包没有伪造 PAT 数据或悄悄替代缺失数据；提供可复现的去模糊与小规模流式CT，并以规模试验考察资源特点。CT是本包的像素驱动小规模离散化，与AIRToolsII不是同一离散矩阵。这些是**参照原文思想的实验**，不是承诺复现原文每个表格数值。

## 3 七种方法的含义与公平性

- **FULL**：同一模型、同一初始化、同一新梯度的无压缩 MM-GKS。允许空间增长；工作缓存按需分段增加。它是最直接的同目标基线。
- **CAPPED**：与FULL相同，但达到同一`kmax`后停止，不压缩。它通常只运行较少步，用于“相同存储上限”比较。不要把它与LM的较长运行误读成相同迭代预算。
- **RESTART**：达到`kmax`后只保留当前解与完整梯度，不保留历史方向。其压缩后维数至多2，与LM默认至多`kmin`不同。`kmin=2`的LM没有历史方向，应与它等价，测试覆盖这一点。
- **TSVD/RBD/SOC/SEC**：有限内存LM，分别使用4种历史选择策略，压缩后仍优先保留当前解和完整梯度。

除 CAPPED 因容量提前停止外，方法使用相同总投影预算；达到梯度精度的正常提前停止仍允许。**不强行保证LM的最终RRE低于充分收敛的FULL。**合理的优势是同预算更好的精度，或达到近似精度所需更少资源。实际结果可以与预期不同，代码不会选择性隐藏失败的方法。

本包的FULL没有冒充原文“每步GCV、旧梯度扩展”的原始软件。这里刻意固定相同目标与新梯度以隔离压缩影响。`robustness`内的(2,1)是固定参数模型基线，不是原文完整GCV复现。

## 4 初始值和算法对应

1. 所有默认方法使用相同`x0=zeros(n,1)`。
2. 从初始残差计算`wd=(e.^2+epsP^2).^(p/2-1)`、`wr=(u.^2+epsQ^2).^(q/2-1)`；实际用`hypot`提高稳定性。
3. 对`D0*A,D0*d`做GKB，将初始解和梯度加入种子空间。
4. 首次投影求出一个近似解，更新双权重。
5. 以`H1=A'*D1^2*A+lambda*Psi'*P1^2*Psi`、`b1=A'*D1^2*d`构造正式Krylov种子。二者均不显式形成。
6. 后续Enlarge更新双权重和完整梯度；达到容量时Compress。

首次投影解**计入**`maxIter`，所有初始化算子调用与时间也计入求解器统计。原文的扩展计数不包含一次性初始化，因此不要混用横坐标口径。

## 5 参数选择是先选后固定

默认`lambdaPolicy='validation'`。在合成实验中，程序生成同一物体的第二次独立含噪观测：

- 每个候选lambda先用训练观测做一个短的LM-TSVD试算；
- 用预测观测与独立验证观测之间的误差评分；Gaussian使用MSE，其余使用平均绝对误差；
- 不使用真实图像或RRE选参数；
- 选择后，FULL与4种LM策略等全部使用同一个lambda，主迭代中不再更新。

此规则假设有独立重复观测；现实只有一次观测时不能声称自动具备验证数据。可以显式采用固定参数：

```matlab
c = struct('lambdaPolicy','fixed','lambdaFixed',0.03);
run_lplq_experiments('compression','quick',c);
```

`lambdaFixed=0.03`只是试算起点，不是保证最佳效果。默认候选网格可覆盖多种量级；当最优点在网格端点时会打印提示，检查`lambda_grid.csv`再扩展网格。验证试算的时间单独记录在`lambda_calibration.mat`，不混入每个正式求解器的计时行；各模型若分别调参，研究整体耗时还应加上它们的调参成本。

对`p!=1`，候选网格还按实验中已知噪声RMS的`p-1`次幂缩放，以免直接沿用不同模型下相同量级的lambda。这是明确的网格设计，不是新的理论选参规则。

`epsP,epsQ`默认都为`1e-3`，一次求解保持固定；已归一化图像且差分不含网格因子。需要敏感性实验时显式更改配置，不能把不同epsilon得到的目标值当作同一函数的下降曲线。

## 6 低内存实现与统计口径

- FFT周期卷积：不生成`n×n`模糊矩阵。
- 周期差分及精确伴随：不生成大规模梯度矩阵。
- 不生成权重对角阵、高维Hessian、高维QR的Q。
- 分块QR处理`[D*AV,D*d]`；用R因子的最后一列保留加权右端。它与正文使用`Qd'*D*d`的最小二乘问题等价，可能多保留一个只影响残差常数的行。
- 只对新增方向调用`A`和`Psi`；压缩优先复用`AV,PsiV`的线性组合。
- 默认只保存标量迭代历史与每种方法的最终解，不保存`saveX(:,iteration)`。
- 每次运行建立唯一结果目录，绝不覆盖之前结果。

**三个不同的内存量：**

1. `activeCacheMiB`：当前有效V/AV/PsiV列所占字节；
2. `allocatedCacheMiB`与`peakAllocatedCacheMiB`：实际分配的三份缓存数组大小，包含预留空列；
3. `estimatedWorkMiB`：含增长/压缩副本、分块QR及向量的工作数组估算。

这些量都不是操作系统测得的 MATLAB 进程峰值。代码没有把缓存估算标成真实进程内存。`memoryBudgetMiB=256`在大型缓存分配前做估算拒绝；这不是操作系统内存隔离，实际进程还包含MATLAB本身、库和用户其他工作区。因此先使用默认quick；不要直接把`imageSize`改成500或1024。内置图像尺寸超过256需要显式放开尺寸开关；CT演示另有96像素/60角度限制。

## 7 HST图像和更大的安全测试

默认明确使用内置phantom，便于无依赖运行。若使用你的HST图像：

```matlab
c = struct();
c.imageSource = 'hst';
c.hstPath = '/Users/zhouting/Documents/编程/文献代码库/L-MM-GKS/LM-MM-GKS/utilities/HSTgray.jpg';
c.imageSize = 96;%(可更改)
%c.p=1;c.q=1;
run_lplq_experiments('compression','quick',c);
```

也可以将`HSTgray.jpg`放入本包`data/`。文件不存在会报错，不悄悄替换图像。HST归一化到[0,1]并缩放，真实数据与`A*xTrue`使用同一个周期边界算子，避免裁边后算子与数据不一致。

较长预算的小规模测试：

```matlab
run_lplq_experiments('compression','paper_small');
```

该profile为96×96、150次求解、`kmax=25`，仍是低内存实验，不是原文500×500、300次的逐项复现。

## 8 输出文件如何看

- `summary.csv`：每个方法的RRE、PSNR、SSIM/HP、实际迭代、状态、时间、缓存、工作数组估算、算子次数等。
- `results.mat`：最终解、所有标量历史、参数及实验数据，方便重新画图。
- `curves.png`：RRE对次数、RRE对时间、分配缓存对次数、完整梯度对次数。
- `memory/scaling.png`：规模实验中随未知量数变化的RRE、计时、缓存和工作数组估算。
- `objective.png`：同一固定目标函数的下降曲线。
- `reconstructions.png`、`absolute_errors.png`：统一色轴重建与绝对误差图。
- `lambda_grid.csv`、`lambda_calibration.mat`：验证选参分数、选定参数和试算总耗时。
- `*_final.png`：仅用于显示的最终图；不会反向影响求解或RRE。

RRE、PSNR使用未截断重建。显示色轴是[0,1]，并未对主迭代施加非负约束。若有兼容外部`HaarPSI`，使用同尺寸二维图像和[0,255]值域；仅该感知评价的重建输入做显示裁剪，并在来源字段标明。缺少HaarPSI或SSIM时填NaN，不用其他指标冒名替代。可手动把你已有的指标工具箱加入MATLAB路径。

`max_iter`表示预算耗尽，不代表梯度已经达到容差；`memory_limit`表示CAPPED容量停止；`converged_gradient`才表示最终真实算子检查满足一阶阈值。小步长仅作为诊断，不被误标为收敛。

## 9 文件用途

完整清单见[docs/FILE_MAP.md](docs/FILE_MAP.md)，参数来源见[docs/PARAMETERS.md](docs/PARAMETERS.md)，实验解读见[docs/EXPERIMENT_GUIDE.md](docs/EXPERIMENT_GUIDE.md)，验证结果见[docs/VALIDATION.md](docs/VALIDATION.md)，模型和论文差异见[docs/ALGORITHM_NOTES.md](docs/ALGORITHM_NOTES.md)。

`examples/`附带实际运行的汇总表、部分图和小规模数据，不需要先运行所有实验才能查看。它们是GNU Octave 8.4.0上的结果；本次没有可用的MATLAB许可证，未声称已在MATLAB或你的Mac上执行。

直接调用求解器示例：

```matlab
setup_lplq;
cfg = lplq_experiment_defaults('quick');
prob = lplq_make_deblur(cfg);
o = lplq_solver_options(cfg,0.03,prob.xTrue);
[xLM,infoLM] = lplq_lmmgks(prob.op,prob.d,o);
[xFull,infoFull] = lplq_mmgks(prob.op,prob.d,o);
```

固定正参数下，压缩理论要求保留当前解和完整梯度；代码会检查保留误差，并以双精度进行核心计算。本包不宣称浮点实现自动等同于精确算术证明，也不保证非凸情形的全局最优或算法一定降低RRE。
