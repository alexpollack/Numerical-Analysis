RGB = imread('pup.tiff');
Isave=imread('pup.tiff');
I1 = rgb2gray(RGB);
I1D=double(I1);

%Compute SVD
[U,S,V]=svd(I1D);
%Finds the original number of singular values
s = diag(S);
rankA = nnz(s);
rankA
for N=5:150:rankA
    %Store the singular values
    C = S;
    %Removes the values not required
    C(N+1:end,:)=0;
    C(:,N+1:end)=0;
    %Image with decreased singular values
    D=U*C*V';
    figure;
    buffer = sprintf('Image with %d singular values', N);
    imshow(uint8(D));
    title(buffer);
end
%Original Image
figure;
imshow(RGB);
title(sprintf('Original image with %f singular values',rankA));