function  fa3D = computeFAfile(filename )
[~,~,raw] = xlsread(filename);

option.compute_eigvec = 0;
for i=2%:size(raw,1)
    
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
        
        S = cat(4,S,S);
        [nx, ny, nz, nu ] = size(S);
        L=nx*ny*nz;
        S_double = reshape(S, [L,nu]);
        b=[b;b];
        [fa,All_eigs]=computeFA(u, b, S_double,option);
        fa3D = reshape(fa, [nx, ny, nz]);
        
       
    end
end
