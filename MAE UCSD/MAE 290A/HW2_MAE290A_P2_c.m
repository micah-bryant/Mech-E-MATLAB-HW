f = F3(:) - mAVec;
fApprox = (U*U'*f) + mAVec;
F = reshape(fApprox,64,64);

subplot(1,2,1)
imagesc(F3); colormap(gray)
title('Original')
subplot(1,2,2)
imagesc(F); colormap(gray)
title('Reconstructed')