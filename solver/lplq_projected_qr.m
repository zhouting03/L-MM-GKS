function [z,fac] = lplq_projected_qr(S,st,d,o)
% 分块 QR 压缩 [D*AV,D*d]：保留加权右端，不形成高维 Q 或正规方程。
k=S.k; Rdaug=zeros(0,k+1); Rr=zeros(0,k);
for lo=1:o.qrBlockRows:numel(d)
    ii=lo:min(numel(d),lo+o.qrBlockRows-1);
    B=[bsxfun(@times,S.Y(ii,1:k),sqrt(st.wd(ii))), sqrt(st.wd(ii)).*d(ii)];
    Rdaug=small_r([Rdaug;B]);
end
for lo=1:o.qrBlockRows:numel(st.u)
    ii=lo:min(numel(st.u),lo+o.qrBlockRows-1);
    B=bsxfun(@times,S.Z(ii,1:k),sqrt(st.wr(ii)));
    Rr=small_r([Rr;B]);
end
fac.Rd=Rdaug(:,1:k); fac.cd=Rdaug(:,k+1); fac.Rr=Rr;
fac.H=[fac.Rd;sqrt(o.lambda)*Rr];
rhs=[fac.cd;zeros(size(Rr,1),1)];
% 小矩阵 SVD 稳定处理可能的数值秩亏；绝不求逆高维 Hessian。
[U,D,V]=svd(fac.H,'econ'); sigma=diag(D);
threshold=max(size(fac.H))*eps*max([sigma;realmin]);
keep=sigma>threshold; fac.rankDeficient=sum(keep)<k;
if any(keep), z=V(:,keep)*((U(:,keep)'*rhs)./sigma(keep));
else, z=zeros(k,1); end
end
function R=small_r(B)
% MATLAB 单输出为 R；Octave 另把反射向量放在下三角，显式去掉它们。
R=triu(qr(B,0)); R=R(1:min(size(R)),:);
end
