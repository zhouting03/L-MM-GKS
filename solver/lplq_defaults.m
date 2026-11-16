function o = lplq_defaults(opts,n)
% 固定模型参数；所有选项显式验证，避免拼写错误被静默忽略。
o.p=1; o.q=1; o.lambda=0.03; o.epsP=1e-3; o.epsQ=1e-3;
o.mode='LM'; o.compression='TSVD'; o.kmin=5; o.kmax=15;
o.maxIter=60; o.init='two-stage'; o.initDim=15; o.x0=[]; o.V0=[];
o.xTrue=[]; o.gradRelTol=1e-6; o.gradAbsTol=1e-10;
o.orthTol=1e-12; o.cacheTol=1e-10; o.descentTol=5e-12;
o.qrBlockRows=1024; o.memoryBudgetMiB=256; o.refreshEvery=10;
o.checkAdjoint=true; o.returnBasis=false; o.verbose=false;
o.rbdTol=1e-5; o.rbdSeed=17; o.socTol=0;
o.secRho=[]; o.secEps=1e-4; o.secMaxIter=12; o.secTol=1e-5;
if nargin<1 || isempty(opts), opts=struct(); end
names=fieldnames(opts);
for j=1:numel(names)
    if ~isfield(o,names{j}), error('LpLq:Option','Unknown option: %s',names{j}); end
    o.(names{j})=opts.(names{j});
end
o.mode=upper(char(o.mode)); o.compression=upper(char(o.compression));
if ~any(strcmp(o.mode,{'LM','FULL','CAPPED','RESTART'}))
    error('LpLq:Option','mode must be LM, FULL, CAPPED or RESTART.');
end
if ~any(strcmp(o.compression,{'TSVD','RBD','SOC','SEC'}))
    error('LpLq:Option','Unknown compression.');
end
if ~any(strcmp(o.init,{'two-stage','simple'})), error('LpLq:Option','Unknown init.'); end
for c={'p','q'}
    validateattributes(o.(c{1}),{'numeric'},{'real','finite','scalar','>',0,'<=',2});
end
for c={'lambda','epsP','epsQ','orthTol','cacheTol','memoryBudgetMiB','secEps'}
    validateattributes(o.(c{1}),{'numeric'},{'real','finite','scalar','positive'});
end
for c={'gradRelTol','gradAbsTol','descentTol','socTol','rbdTol','secTol'}
    validateattributes(o.(c{1}),{'numeric'},{'real','finite','scalar','nonnegative'});
end
for c={'kmin','kmax','maxIter','qrBlockRows','refreshEvery','secMaxIter'}
    validateattributes(o.(c{1}),{'numeric'},{'scalar','integer','positive'});
end
validateattributes(o.initDim,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(o.rbdSeed,{'numeric'},{'scalar','integer','nonnegative'});
if o.kmin<2 || o.kmax<=o.kmin || o.kmax>n
    error('LpLq:Option','Require 2 <= kmin < kmax <= n.');
end
if isempty(o.x0), o.x0=zeros(n,1); end
o.x0=double(o.x0(:));
if numel(o.x0)~=n || ~isreal(o.x0) || any(~isfinite(o.x0))
    error('LpLq:Option','x0 must be a finite real n-vector.');
end
if ~isempty(o.xTrue)
    o.xTrue=double(o.xTrue(:));
    if numel(o.xTrue)~=n || ~isreal(o.xTrue) || any(~isfinite(o.xTrue))
        error('LpLq:Option','xTrue has incorrect size or values.');
    end
end
if ~isempty(o.V0) && (size(o.V0,1)~=n || ~isreal(o.V0) || any(~isfinite(o.V0(:))))
    error('LpLq:Option','V0 has incorrect size or values.');
end
if isempty(o.secRho), o.secRho=o.lambda; end
validateattributes(o.secRho,{'numeric'},{'real','finite','scalar','positive'});
end
