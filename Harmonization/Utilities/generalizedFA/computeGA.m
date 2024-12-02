function  ga3D = computeGA(filename )
[~,~,raw] = xlsread(filename);

for i=1%:size(raw,1)
    
    if  exist(raw{i,1},'file')
        [S, u, b, voxel, Sm, s0, ~, ~] = loadDWI(raw{i,1}, ...
            'nii.gz', 1);
        
        s0(isnan(s0))=0;
        s0(s0<=0)=1;
        
        
        S(S<0+eps) = 0; S(S>1-eps) = 1;
        S(isnan(S)) = 0;
        
        mpath = raw{i,2};
        ind2 = strfind(mpath,'.');
        masktype = mpath(ind2+1:end);
        [mask, sd_mask, ~] = loadMask(raw{i,2}, 'nii.gz', 1);
        S = S.*single(repmat(mask,[1,1,1,size(S,4)]));
        
        
        ga3D = generalizedAnisotropy(S);
    end
end
