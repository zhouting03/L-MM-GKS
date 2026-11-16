# 文件用途对照

所有MATLAB文件都有`%`简要注释。主函数名使用`lplq_`前缀，不覆盖你的原始MMGKS/LMMGKS。

| 文件 | 主要用途 |
|---|---|
| `setup_lplq.m` | 添加本包目录 |
| `run_lplq_experiments.m` | 实验统一入口、唯一输出目录、保存配置与调参记录 |
| `solver/lplq_lmmgks.m` | 固定参数求解主循环、实际计数、最终梯度复核 |
| `solver/lplq_mmgks.m` | 同目标无压缩基线包装 |
| `solver/lplq_defaults.m` | 选项默认值及严格验证 |
| `solver/lplq_memory_plan.m` | 在大型缓存分配前估算并检查工作内存预算 |
| `solver/lplq_state.m` | 当前目标、双权重、未投影完整梯度 |
| `solver/lplq_keep_basis.m` | 优先保留解/梯度的双重MGS |
| `solver/lplq_cache.m` | 初始化V/AV/PsiV缓存 |
| `solver/lplq_append.m` | 按需扩容并添加新增方向 |
| `solver/lplq_projected_qr.m` | 分块双加权QR和小空间SVD求解，不构造高维Q |
| `solver/lplq_gkb_seed.m` | 初始数据加权GKB右基 |
| `solver/lplq_krylov_seed.m` | 第一次求解后的双权重正规算子种子 |
| `solver/lplq_initialize.m` | 完整两阶段初始化或简化初始化 |
| `solver/lplq_step.m` | 一次代理投影及数值下降保障 |
| `solver/lplq_enlarge.m` | 更新新权重和完整梯度并扩展空间 |
| `solver/lplq_select_basis.m` | 真正的tSVD/RBD/SOC/SEC历史选择器 |
| `solver/lplq_compress.m` | 直接保留高维解/梯度、检查保留误差、复用缓存 |
| `solver/lplq_record.m` | 只生成标量诊断历史 |
| `operators/lplq_from_matrices.m` | 包装矩阵或Aclass，避免size(A,2)兼容问题 |
| `operators/lplq_grad2.m` | 周期二维差分与其伴随 |
| `operators/lplq_blur_fft.m` | 矩阵自由周期卷积与伴随 |
| `operators/lplq_motion_psf.m` | 无图像工具箱的运动模糊核 |
| `operators/lplq_phantom.m` | 可复现的内置测试图像 |
| `operators/lplq_noise.m` | Gaussian/Laplace/混合脉冲噪声，保留随机数状态 |
| `operators/lplq_make_deblur.m` | 创建同一正向模型的数据、验证观测与图像 |
| `operators/lplq_tomo_matrix.m` | 带尺寸限制的小规模稀疏CT矩阵 |
| `experiments/lplq_experiment_defaults.m` | quick与paper_small配置 |
| `experiments/lplq_solver_options.m` | 从实验配置建立相同的求解器设置 |
| `experiments/lplq_select_lambda.m` | 单独的验证选参阶段，正式迭代不再改变参数 |
| `experiments/lplq_compare.m` | 逐方法求解、统计、打印与保存 |
| `experiments/lplq_metrics.m` | RRE/PSNR，以及可选真实SSIM/HaarPSI |
| `experiments/lplq_write_csv.m` | 不依赖table的标准CSV写出 |
| `experiments/lplq_plot_results.m` | 同目标对比图和误差图 |
| `experiments/lplq_plot_scale.m` | 跨尺寸RRE、计时、缓存和工作数组估算汇总图 |
| `experiments/lplq_streaming_demo.m` | 分块CT的独立、热启动、回收和全数据参考 |
| `tests/test_lplq_all.m` | 梯度、加权QR、35组模式/模型、退化情况与资源预检 |

建议阅读顺序：README → lplq_lmmgks → lplq_initialize → lplq_enlarge/lplq_step → lplq_compress → lplq_select_basis。
