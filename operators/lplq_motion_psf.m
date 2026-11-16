function P = lplq_motion_psf(width,len,angleDeg)
% 亚像素直线采样构造运动模糊 PSF，不依赖 fspecial/imresize。
if mod(width,2)~=1 || width<3, error('LpLq:PSF','Use an odd PSF width >= 3.'); end
P=zeros(width); center=(width+1)/2;
for t=linspace(-len/2,len/2,401)
    r=center+t*sind(angleDeg); c=center+t*cosd(angleDeg);
    r0=floor(r); c0=floor(c); fr=r-r0; fc=c-c0;
    for a=0:1
        for b=0:1
            rr=r0+a; cc=c0+b;
            if rr>=1 && rr<=width && cc>=1 && cc<=width
                wr=(1-fr)*(a==0)+fr*(a==1); wc=(1-fc)*(b==0)+fc*(b==1);
                P(rr,cc)=P(rr,cc)+wr*wc;
            end
        end
    end
end
P=P/sum(P(:));
end
