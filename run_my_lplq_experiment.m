% 参数集中入口：放入原 LpLq_MM_GKS_MATLAB 根目录，点击“运行”。
% 只需修改下面的普通变量；后面的结构体映射不必手动调整。
% 本入口运行 compression：各方法采用同一 p,q、数据和固定 lambda。

%% 1. 图像与观测：HST 路径按你的电脑设置
imageSource = 'hst';              % 'hst' 或 'phantom'（内置测试图，无需文件）
hstPath = '/Users/zhouting/Documents/编程/文献代码库/L-MM-GKS/LM-MM-GKS/utilities/HSTgray.jpg';
imageSize = 128;                    % 缩放为 imageSize × imageSize；先从64或96开始
noiseType = 'laplace_impulse';     % 'gaussian'、'laplace'、'laplace_impulse'
noiseLevel = 0.02;                 % 基础噪声标准差/干净观测RMS；不含附加脉冲
impulseFraction = 0.02;            % 每个观测点被叠加脉冲的概率，仅混合噪声使用
seed = 17;                        % 固定种子便于不同方法使用同一数据
% 现有代码的脉冲幅度固定为0.5*max(干净观测极差,干净观测RMS)。
% 模糊核固定为9×9、长度7、角度35度；源代码位置见配套说明。

%% 2. 主目标参数：一次正式求解中全部固定
p = 2;                            % 数据项指数；例如改为0.6
q = 1;                            % 正则项指数；例如改为1.2
epsP = 1e-3;                      % 数据项平滑参数
epsQ = 1e-3;                      % 正则项平滑参数

%% 3. lambda：'validation'先验证后固定；'fixed'直接使用lambdaFixed
lambdaPolicy = 'validation';
lambdaFixed = 0.03;               % 只有lambdaPolicy='fixed'时生效
lambdaGrid = [1e-3,3e-3,1e-2,3e-2,1e-1,3e-1];
pilotSteps = 25;                   % 每个候选的试算预算，不计入正式maxIter
% p不等于1时实际网格会乘 max(噪声实际RMS,1e-3)^(p-1)。
% 最优点在最小端时，可将网格改为[1e-4,3e-4,lambdaGrid]后重新验证。

%% 4. 子空间、总预算与主要工作数组限制
kmin = 5;                         % 压缩后至多保留的维数，包含解与完整梯度
kmax = 25;                        % 有限内存方法的空间容量；FULL允许继续增长
initDim = 13;                     % 初始GKB步数，实际截为min(initDim,kmax-2)
maxIter = 150;                     % 全过程投影求解预算，包含首次初始化求解
memoryBudgetMiB = 256;            % 工作数组估算的拒绝阈值，不是进程内存上限
qrBlockRows = 1024;               % 分块QR的行数
methods = {'FULL','CAPPED','RESTART','TSVD','RBD','SOC','SEC'};
% 例如只比较三种：methods = {'FULL','CAPPED','TSVD'};

%% 5. 输出：不再生成FULL_final.png等七张独立小图
makeFigures = true;               % 保留组合重建图、误差图、curves.png等
saveFinalImages = false;          % 关闭每种算法的独立*_final.png

%% 以下自动传递参数；通常不必修改
entryDir = fileparts(mfilename('fullpath'));
if exist(fullfile(entryDir,'setup_lplq.m'),'file') == 2
    addpath(entryDir);
elseif exist('setup_lplq','file') ~= 2
    error('LpLq:EntryLocation', ...
        '请把本脚本放入原代码包根目录，与setup_lplq.m放在一起。');
end
setup_lplq;
c = struct();
c.imageSource=imageSource; c.hstPath=hstPath; c.imageSize=imageSize;
c.noiseType=noiseType; c.noiseLevel=noiseLevel;
c.impulseFraction=impulseFraction; c.seed=seed;
c.p=p; c.q=q; c.epsP=epsP; c.epsQ=epsQ;
c.lambdaPolicy=lambdaPolicy; c.lambdaFixed=lambdaFixed;
c.lambdaGrid=lambdaGrid; c.pilotSteps=pilotSteps;
c.kmin=kmin; c.kmax=kmax; c.initDim=initDim; c.maxIter=maxIter;
c.memoryBudgetMiB=memoryBudgetMiB; c.qrBlockRows=qrBlockRows;
c.methods=methods; c.makeFigures=makeFigures; c.saveFinalImages=saveFinalImages;
fprintf('\n实际请求：p=%g, q=%g, epsP=%g, epsQ=%g\n',p,q,epsP,epsQ);
fprintf('图像：%s，%d×%d；噪声：%s；noiseLevel=%g；seed=%d\n', ...
    imageSource,imageSize,imageSize,noiseType,noiseLevel,seed);
fprintf('kmin=%d, kmax=%d, 总求解预算=%d；lambda策略=%s\n', ...
    kmin,kmax,maxIter,lambdaPolicy);
outdir = run_lplq_experiments('compression','quick',c);

%% 用求解器输出复核参数，避免只修改了未生效的默认值
savedRun = load(fullfile(outdir,'compression','results.mat'),'results','noise','lambda');
for j=1:numel(savedRun.results)
    actual=savedRun.results{j}.info.options;
    assert(actual.p==p && actual.q==q, 'LpLq:ParameterMismatch', ...
        '求解器实际p,q与入口设置不一致，请检查MATLAB路径。');
end
fprintf('\n已复核全部方法：实际p=%g，q=%g，正式lambda=%.12g\n',p,q,savedRun.lambda);
fprintf('实际总噪声相对L2=%.5g（包含脉冲时不等于noiseLevel）\n', ...
    savedRun.noise.actualRelativeL2);
fprintf('本次目录：%s\n',outdir);
if strcmp(lambdaPolicy,'validation')
    calibration = load(fullfile(outdir,'compression','lambda_calibration.mat'),'cal');
    if calibration.cal.onBoundary
        fprintf('提示：验证最优点位于候选端点，请检查lambda_grid.csv并扩展网格。\n');
    end
end
% 关闭saveFinalImages只影响本次新输出，已有目录中的旧PNG不会被删除。
