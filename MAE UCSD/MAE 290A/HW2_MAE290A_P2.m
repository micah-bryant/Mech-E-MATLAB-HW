mAvg = 0; mVec = zeros(4096, 132);

for i = 1:length(Y)
    mTemp = Y{i};
    mAvg = mAvg+mTemp;
    mVec(:,i) = mTemp(:);
end
mAvg = mAvg/length(Y);
mAVec = mAvg(:);
D = (1/sqrt(132)) .* (mVec-mAVec);
[U,S,V] = svd(D,'econ');

subplot(2,2,1)
imagesc(reshape(U(:,1),64,64)); colormap(gray)
title('First Singular Vector')
subplot(2,2,2)
imagesc(reshape(U(:,2),64,64)); colormap(gray)
title('Second Singular Vector')
subplot(2,2,3)
imagesc(reshape(U(:,3),64,64)); colormap(gray)
title('Third Singular Vector')
subplot(2,2,4)
imagesc(reshape(U(:,4),64,64)); colormap(gray)
title('Fourth Singular Vector')