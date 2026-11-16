function prob = lplq_make_deblur(cfg)
% 原文去模糊测试的低内存版本；phantom 默认，HST 必须显式指定。
n=cfg.imageSize;
validateattributes(n,{'numeric'},{'scalar','integer','>=',16});
if n>256 && ~cfg.allowLargerImages
    error('LpLq:SizeGuard','imageSize > 256 requires allowLargerImages=true.');
end
switch lower(cfg.imageSource)
    case 'phantom'
        X=lplq_phantom(n); source='built-in piecewise-constant phantom';
    case 'hst'
        candidates={cfg.hstPath,fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','HSTgray.jpg')};
        path='';
        for j=1:numel(candidates)
            if exist(candidates{j},'file')==2, path=candidates{j}; break; end
        end
        if isempty(path), error('LpLq:HST','HSTgray.jpg not found: set cfg.hstPath. No silent image substitution.'); end
        X=double(imread(path)); if ndims(X)==3, X=mean(X,3); end
        [a,b]=size(X); [xx,yy]=meshgrid(linspace(1,b,n),linspace(1,a,n));
        X=interp2(X,xx,yy,'linear'); X=X-min(X(:)); X=X/max(max(X(:)),realmin); source=path;
    otherwise
        error('LpLq:Image','imageSource must be phantom or hst.');
end
P=lplq_motion_psf(9,7,35); padded=zeros(n); padded(1:size(P,1),1:size(P,2))=P;
H=fft2(circshift(padded,-floor(size(P)/2)));
op.n=n*n; op.m=n*n; op.r=2*n*n;
op.A=@(v) lplq_blur_fft(v,H,n,n); op.At=@(v) lplq_blur_fft(v,conj(H),n,n);
op.Psi=@(v) lplq_grad2(v,n,n,false); op.Psit=@(v) lplq_grad2(v,n,n,true);
truth=X(:); clean=op.A(truth);
[d,meta]=lplq_noise(clean,cfg.noiseType,cfg.noiseLevel,cfg.seed,cfg.impulseFraction);
[dval,~]=lplq_noise(clean,cfg.noiseType,cfg.noiseLevel,cfg.seed+1009,cfg.impulseFraction);
prob.op=op; prob.d=d; prob.validationData=dval; prob.xTrue=truth;
prob.imageSize=[n,n]; prob.imageSource=source; prob.psf=P; prob.noise=meta;
prob.clean=clean;
end
