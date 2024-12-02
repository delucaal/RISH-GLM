% load DWI data
% DWI can be in three formats: {'nhdr'}, {'mat'} or
% {'nii','img','nii.gz','img.gz'}

%% Inputs:
% imagefullpath : full path DWI data file
% imagetype: image type of DDWI data file

%% Outputs:
 % S - the DWI data in 4D
 % u - the gradient directions
 % b - the b-value
 % voxel - the voxel size
 % Sm - the space directions matrix
function [S, u, b, voxel, Sm, s0, s0Arr, option] = loadDWI(imagefullpath, imagetype)
    % This part has been completely rewritten by Alberto De Luca
    
    fprefix = imagefullpath(1:end-length(imagetype));
    S = load_untouch_nii([fprefix imagetype]);
    option = S;
    S = single(S.img);%*S.hdr.dime.scl_slope+S.hdr.dime.scl_inter;
    b = load([fprefix 'bval'])';
    u = load([fprefix 'bvec'])';

    b0s = b<=1;
    b(b0s) = 0;
    s0 = mean(S(:,:,:,b0s),4);
    S = S(:,:,:,~b0s);
    for ij=1:size(S,4)
       V = S(:,:,:,ij);
       V = V./s0;
       V(s0 == 0) = 0;
       S(:,:,:,ij) = V;
    end

    b = b(~b0s);
    u = u(~b0s,:);
    u = [u; -u];

    voxel = option.hdr.dime.pixdim(2:4);
    Sm = [];
    s0Arr = [];
    
    return
    %nnrd
    if strcmp(imagetype, 'nhdr') || strcmp(imagetype, 'nrrd') 
        
        dwi = nrrdLoadWithMetadata(imagefullpath); 
        S =  dwi.data;
        
        u =  dwi.gradientdirections;
        permuted = 1;
        switch size(u,1)
            case size(S,4)
               order = [1 2 3 4];
               permuted =0;
            case size(S,1)
                order = [2 3 4 1];
            case size(S,2)
                order = [1 4 3 2];
            case size(S,3)
                order = [1 2 4 3];
        end
        S = permute(S,order);
        dwi.data = S;
        voxel = abs(diag(dwi.spacedirections));
        
        u_sum = sqrt(sum(u.^2,2));
        ind_0 =  u_sum <= 0.05;
        u(ind_0, :) = [];
        u_sum(ind_0,:)=[];
        b = repmat(dwi.bvalue, length(u), 1);
        
        
        s0 = mean(S(:,:,:,ind_0),4);       
        s0Arr = S(:,:,:,ind_0);
        S(:,:,:,ind_0)=[];
        s0(s0==0) = 1;
        
        S = single(S);
        Sm=[dwi.spacedirections;dwi.spaceorigin'];
        
        for j=1:size(S,4)
           S(:,:,:,j) = S(:,:,:,j)./(s0+eps);
        end
        
      
        
        u=u./repmat(u_sum,1, 3);
        
        
        u = [u; -u]; 
        
        option = dwi;
        
                
    %nifti 
    elseif(strcmp(imagetype, 'nii') || strcmp(imagetype, 'nii.gz')...
            || strcmp(imagetype, 'img') || strcmp(imagetype, 'img.gz')) %single shell only
        nii = load_untouch_nii(imagefullpath); voxel = nii.hdr.dime.pixdim(2:4);
      
        Sm=[nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z]';
        S = nii.img; option= nii; clear nii;    
        
        
        upath = [imagefullpath(1:end-length(imagetype)) 'bvec'];
        fid = fopen(upath,'r'); u = fscanf(fid, '%f',[inf ]); fclose(fid);
        
        %ARISTOTLE data 3xnumber of gradients
        [status,output]=system(['wc -l < ', upath]);
        if status==0 & str2num(output)==3
            usz=length(u)/3;
            uhat=zeros(usz,3);
            for i=1:usz
                uhat(i,1)=u(i);
                uhat(i,2)=u(i+usz);
                uhat(i,3)=u(i+2*usz);
            end
            u=uhat;
        else
            usz=length(u)/3;
            uhat=zeros(usz,3);
            ccc=1;
            for i=1:usz
                uhat(i,1)=u(ccc,1);
                uhat(i,2)=u(ccc+1,1);
                uhat(i,3)=u(ccc+2,1);
                ccc=ccc+3;
            end
            u=uhat;
        end
        
        bpath = [imagefullpath(1:end-length(imagetype)) 'bval'];
        fid = fopen(bpath,'r'); b = fscanf(fid, '%f'); fclose(fid);   
        
%         %%FOR CDMRI CHALLENGE: I removed first b0, correct it later for the
%         %%other data
%         S(:,:,:,1)=[];
%         u(1,:)=[];
%         b(1)=[];
%         %%
%        
        u_sum = sqrt(sum(u.^2,2));
       
        
        ind_0 =  u_sum <= 0.05 | b==0;
        s0 = mean(S(:,:,:,ind_0),4);
       
        s0Arr = S(:,:,:,ind_0);
        S(:,:,:,ind_0)=[];
        
        
        S = single(S);
       
        %% FOR CDMRI CHALLENGE:I removed all b0s to observe the effects of the non normalized data
        %maxs0 = max(s0(:));
        %S(S>maxs0)=maxs0;
        for j=1:size(S,4)
           S(:,:,:,j) = S(:,:,:,j)./(s0+eps);
        end
        %s0(:)=1;
       %%
         
        u(ind_0,:)=[];
        u_sum(ind_0,:)=[];
        u=u./repmat(u_sum,1, 3);
        u = [u; -u]; 
        b(ind_0)=[];
        
       
       
    else
        fprintf(fid,'Wrong file type! Image should be in .nrrd or in .mat type');
        S = []; u=[]; b=[]; voxel = []; s0 = []; Sm=[];
    end

end