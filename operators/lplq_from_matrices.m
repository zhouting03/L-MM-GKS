function op = lplq_from_matrices(A,L)
% 显式小矩阵或支持 A*x、A'*x 的 Aclass 包装；不用 size(A,2)。
[m,n]=size(A); [r,nL]=size(L);
if nL~=n, error('LpLq:Dimension','A and L have different column counts.'); end
op.n=n; op.m=m; op.r=r;
op.A=@(x) A*x; op.At=@(y) A'*y;
op.Psi=@(x) L*x; op.Psit=@(y) L'*y;
end
