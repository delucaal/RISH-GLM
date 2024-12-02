function  ga3D = computeGeneralizedAnisotropy(dwifile, dwimask, outputGFAfile)


    if  exist(dwifile,'file')	
	% read dwifile /bvec and bval is assumed to share the same name with dwi file e.g if the file name is  dwi.nii.gz, then  dwi.bval dwi.bvec should be in the same folder.
        [S, u, b, voxel, Sm, s0, ~, ~] = loadDWI(dwifile, ...
            'nii.gz', 1);
        
        s0(isnan(s0))=eps;
        s0(s0<=0)=eps;
        
        
        S(S<0+eps) = eps; S(S>1-eps) = 1;
        S(isnan(S)) = eps;
        
       
       % read mask.
        [mask, sd_mask, nii] = loadMask(dwimask, 'nii.gz', 1);
        S = S.*single(repmat(mask,[1,1,1,size(S,4)]));
        
	% GFA        
        ga3D = generalizedAnisotropy(S);

	% save GFA
	nii.img = ga3D;
	nii.hdr.dime.datatype=16;
	save_untouch_nii(nii, outputGFAfile);
    end
end
