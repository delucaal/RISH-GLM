function rishFeatures_file(dwifile, maskfile, outputPrefix)
addpath(genpath('./Utilities'));
SHorder= 8; lambda=0.006;
%load dwi data
[S, u, b, voxel, Sm, s0, ~, option] = loadDWI(dwifile, ...
   'nii.gz');

s0(isnan(s0))=0;
s0(s0<=0)=1;


%Remove possible noise
S(S<0+eps) = 0;
S(S>1-eps) = 1;
S(isnan(S)) = 0;


%load brain mask
[mask, sd_mask, moption] = loadMask(maskfile, 'nii.gz', 1);
S = S.*single(repmat(mask,[1,1,1,size(S,4)]));
u = u(1:size(S,4),:);

% compute rish features
IMG = rishFeatures(S, u,  mask, SHorder, lambda);          


for l=1:SHorder/2+1
    moption.img = squeeze(IMG.rishImgs(:,:,:,l));
    moption.filetype = 16;
    moption.hdr.dime.datatype = 16;
   

    save_untouch_nii(moption, [outputPrefix, '_L', num2str((l-1)*2), '.nii.gz']);

end
end