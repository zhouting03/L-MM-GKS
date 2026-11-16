function y = lplq_blur_fft(x,H,nx,ny)
% 周期卷积；伴随调用传入 conj(H)，避免显式 n-by-n 模糊矩阵。
y=real(ifft2(H.*fft2(reshape(x,nx,ny)))); y=y(:);
end
