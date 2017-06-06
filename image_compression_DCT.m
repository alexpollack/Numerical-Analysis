RGB = imread('pup.tiff');
Isave=imread('pup.tiff');
I = rgb2gray(RGB);
I = im2double(I);

T = dctmtx(8);
dct1 = @(block_struct) block_struct.data * T';
Bt1 = blockproc(I,[8 8],dct1);
dct2 = @(block_struct) T * block_struct.data;
Bt2 = blockproc(Bt1,[8 8],dct2);

dct12 = @(block_struct) block_struct.data * T';
Bt12 = blockproc(I,[8 8],dct12);
dct22 = @(block_struct) T * block_struct.data;
Bt22 = blockproc(Bt1,[8 8],dct22);

dct13 = @(block_struct) block_struct.data * T';
Bt13 = blockproc(I,[8 8],dct13);
dct23 = @(block_struct) T * block_struct.data;
Bt23 = blockproc(Bt1,[8 8],dct23);

mask = [1   1   1   0   0   0   0   0;
        1   1   0   0   0   0   0   0;
        1   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0;
        0   0   0   0   0   0   0   0];
mask2 = [1   1   1   1   1   0   0   0;
         1   1   1   1   0   0   0   0;
         1   1   1   0   0   0   0   0;
         1   1   0   0   0   0   0   0;
         1   0   0   0   0   0   0   0;
         0   0   0   0   0   0   0   0;
         0   0   0   0   0   0   0   0;
         0   0   0   0   0   0   0   0];
mask3 = [1   1   1   1   1   1   1   0;
         1   1   1   1   1   1   0   0;
         1   1   1   1   1   0   0   0;
         1   1   1   1   0   0   0   0;
         1   1   1   0   0   0   0   0;
         1   1   0   0   0   0   0   0;
         1   0   0   0   0   0   0   0;
         0   0   0   0   0   0   0   0];
B21 = blockproc(Bt2,[8 8],@(block_struct) mask .* block_struct.data);
B22 = blockproc(Bt22,[8 8],@(block_struct) mask2 .* block_struct.data);
B23 = blockproc(Bt23,[8 8],@(block_struct) mask3 .* block_struct.data);

invdct1 = @(block_struct) block_struct.data * T;
It2 = blockproc(B21,[8 8],invdct1);
invdct2 = @(block_struct) T' * block_struct.data;
I2 = blockproc(It2,[8 8],invdct2);

invdct12 = @(block_struct) block_struct.data * T;
It22 = blockproc(B22,[8 8],invdct12);
invdct22 = @(block_struct) T' * block_struct.data;
I22 = blockproc(It22,[8 8],invdct22);

invdct13 = @(block_struct) block_struct.data * T;
It23 = blockproc(B23,[8 8],invdct13);
invdct23 = @(block_struct) T' * block_struct.data;
I23 = blockproc(It23,[8 8],invdct23);

figure
imshow(Isave)
title(sprintf('Original'));
figure
imshow(I23)
title(sprintf('Mask 3'));
figure
imshow(I22)
title(sprintf('Mask 2'));
figure
imshow(I2)
title(sprintf('Mask 1'));
