function [d,meta] = lplq_noise(clean,kind,level,seed,impulseFraction)
% 按分布尺度加噪；不在采样后强行二范数归一化或裁剪观测。
if nargin<5, impulseFraction=0.02; end
old=rng; rng(seed,'twister');
clean=clean(:); m=numel(clean); rmsSignal=norm(clean)/sqrt(m);
sigma=level*max(rmsSignal,realmin);
kind=lower(char(kind)); noise=zeros(m,1); laplaceScale=0;
switch kind
    case 'gaussian'
        noise=sigma*randn(m,1);
    case {'laplace','laplace_impulse'}
        laplaceScale=sigma/sqrt(2); u=rand(m,1)-0.5;
        noise=-laplaceScale*sign(u).*log(max(1-2*abs(u),realmin));
    otherwise
        rng(old); error('LpLq:Noise','Unknown noise model.');
end
actualFraction=0;
if strcmp(kind,'laplace_impulse')
    mask=rand(m,1)<impulseFraction;
    amplitude=0.5*max(max(clean)-min(clean),rmsSignal);
    noise(mask)=noise(mask)+amplitude*sign(randn(sum(mask),1));
    actualFraction=mean(mask);
end
rng(old); d=clean+noise;
meta.kind=kind; meta.requestedRMSFraction=level;
meta.actualRelativeL2=norm(noise)/max(norm(clean),realmin);
meta.rms=norm(noise)/sqrt(m); meta.laplaceScale=laplaceScale;
meta.impulseFraction=actualFraction; meta.seed=seed;
end
