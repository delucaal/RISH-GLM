function [highResImg, highResMask, sp_high, Sm_n, s0_n,  so_mask, sd_mask]=resampleData(lowResImg, lowResMask, sp_low, sp_high, Sm, s0, sd_mask, splineorder)

[sz_x, sz_y, sz_z, n]=size(lowResImg);
step(1) = sp_high(1)/sp_low(1);
step(2) = sp_high(2)/sp_low(2);
step(3) = sp_high(3)/sp_low(3);


sOrder = splineorder;
%[x,y,z]=ndgrid(1-step(1)/2:step(1):sz_x+step(1)/2, 1-step(2)/2:step(2):sz_y+step(2)/2, 1-step(3)/2:step(3):sz_z+step(3)/2);
[x,y,z]=ndgrid(1:step(1):(sz_x+step(1)+0.01), 1:step(2):(sz_y+step(2)+0.01), 1:step(3):(sz_z+step(3)+0.01));
disp(size(x))
highResImg=zeros(size(x,1) , size(x,2), size(x,3), n);
for i=1:n
   v1=spm_bsplinc(double(squeeze(lowResImg(:,:,:,i))),[sOrder sOrder sOrder 0 0 0]);
   highResImg(:,:,:,i) = single(spm_bsplins(v1,x,y,z,[sOrder sOrder sOrder 0 0 0]));
   i
end
%matlabpool close;
%poolobj = gcp('nocreate');
%delete(poolobj);

ind = find(lowResMask>0);

highResMask = interp3(1:sz_y, 1:sz_x, 1:sz_z,double( lowResMask), y, x, z,'linear');
highResMask = double(highResMask>=0.5);
se = strel3d(3);
highResMask = imopen(highResMask, se);

sOrder =splineorder;
v1 = spm_bsplinc(double(s0),[sOrder sOrder sOrder 0 0 0]);
s0_n = spm_bsplins(v1,x,y,z,[sOrder sOrder sOrder 0 0 0]);
s0_n = s0_n.*highResMask;
s0_n(s0_n<=0) = 1;
s0_n(isnan(s0_n))=1;
deGibbs = RemoveGibbsRinging(s0_n);
s0_n =  deGibbs;

s0_n(s0_n>max(s0(:))) = max(s0(ind));
s0_n(s0_n<min(s0(:))) = min(s0(ind));

s0_n(s0_n<=0) = 1;
s0_n(isnan(s0_n))=1;



lowResImg_r = reshape(lowResImg, sz_x*sz_y*sz_z, n);
highResImg(highResImg>max(lowResImg(:))+eps) = max(max(lowResImg_r(ind,:)));
highResImg(highResImg<min(lowResImg(:))-eps) = min(min(lowResImg_r(ind,:)));

for j=1:size(highResImg,4)
   highResImg(:,:,:,j) = highResImg(:,:,:,j)./s0_n;
   highResImg(:,:,:,j) = highResImg(:,:,:,j).*highResMask;
end


highResImg(highResImg(:)<0+eps) = 0;
highResImg(highResImg(:)>1-eps) = 1;
highResImg(isnan(highResImg)) = 0;







% for dwi
Sm_n = diag([sign(Sm(1,1))*sp_high(1)   sign(Sm(2,2))*sp_high(2) sign(Sm(3,3))*sp_high(3) ]);
sz = size(x);
Sm_n(4,:) = -0.5 * Sm_n * sz(1:3)'; %the space origin in Slicer

% for mask
sd_mask = diag([sign(sd_mask(1,1))*sp_high(1)   sign(sd_mask(2,2))*sp_high(2) sign(sd_mask(3,3))*sp_high(3) ]);
so_mask = -0.5 * sd_mask * sz(1:3)'; %the space origin in Slicer
disp('End of resampling...');
end