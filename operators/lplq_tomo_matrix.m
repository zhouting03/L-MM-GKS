function [A,nDet] = lplq_tomo_matrix(n,angles)
% 小规模平行束的像素驱动离散化；不是 AIRToolsII 的严格相同离散矩阵。
% 每像素每角度最多两个非零贡献，构造前限制规模。
if n>96 || numel(angles)>60, error('LpLq:TomoGuard','Use n<=96 and <=60 angles for the bundled sparse CT demo.'); end
nDet=ceil(sqrt(2)*n)+2; N=n*n; na=numel(angles);
[xx,yy]=meshgrid((1:n)-(n+1)/2,(1:n)-(n+1)/2); xx=xx(:); yy=yy(:);
rows=zeros(2*N*na,1); cols=rows; vals=rows; used=0;
for j=1:na
    t=xx*cosd(angles(j))+yy*sind(angles(j))+(nDet+1)/2;
    lo=floor(t); alpha=t-lo;
    for s=0:1
        bins=lo+s; good=bins>=1 & bins<=nDet;
        cc=find(good); len=numel(cc); ii=used+(1:len);
        rows(ii)=(j-1)*nDet+bins(good); cols(ii)=cc;
        if s==0, vals(ii)=1-alpha(good); else, vals(ii)=alpha(good); end
        used=used+len;
    end
end
A=sparse(rows(1:used),cols(1:used),vals(1:used),nDet*na,N);
end
