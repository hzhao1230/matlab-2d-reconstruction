function image = remove_single(image)
% -------------------------------------------------------------------------
% Eliminate all the single pixels: isolated "1" or "0";
% Used in 2Dconstruction_with_iteration.m,  remove_double.m
% -------------------------------------------------------------------------
[x,y]=size(image);
cz=zeros(x,1);
lz=zeros(1,y);

t=x*y+1;

% (1) remove all the isolated "1" pixels:
img1=[cz, image]; img1=img1(1:x,1:y);
img2=image(1:x,2:y); img2=[img2, cz];
img3=[lz; image]; img3=img3(1:x,1:y);
img4=image(2:x,1:y); img4=[img4; lz];
img=(image-img1)+(image-img2)+(image-img3)+(image-img4);
img_one = abs( ceil( abs( (img-4)/t )  ) );  img_one = abs(img_one-1);
image = image-img_one;

% (2) remove all the isolated "0" pixels:
img1=[cz, image]; img1=img1(1:x,1:y);
img2=image(1:x,2:y); img2=[img2, cz];
img3=[lz; image]; img3=img3(1:x,1:y);
img4=image(2:x,1:y); img4=[img4; lz];
img=(image-img1)+(image-img2)+(image-img3)+(image-img4);
img_zero = abs( ceil( abs( (img+4)/t )  ) );  img_zero = abs(img_zero-1);
image = image+img_zero;