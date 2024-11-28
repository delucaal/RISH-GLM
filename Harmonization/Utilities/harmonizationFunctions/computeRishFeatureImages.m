%This function computes RISH feature images (L0, L2, L4, L6, L8) of each data and saves the images under rishPath

function SITES = computeRishFeatureImages(SITES, coption)

SHorder = coption.SHOrder; % Spherical Harmonics Order
fid1 = fopen(['./log/log_', coption.name,'computerRishFeatureImages'],'w');


for s=1:length(SITES)
    for i=1:SITES{s}.ImageNum
        data =  SITES{s}.InputImages{i};
        fprintf(fid1,'%s\n', data.ImageFullPath);
        
        
        [status, ~] =system(['ls  ', data.Image_Harmonized_dir 'rish/' data.ImageName, '*_L0.nii.gz  > ', data.Image_Harmonized_dir 'data.txt']);
        %if ~exist([data.Image_Harmonized_dir 'rish/' data.ImageName newfile_suffix  '_L0.nii.gz'],'file') || coption.force
        if (status) || coption.force  || ~exist([data.Image_Harmonized_dir,'newfile_suffix.mat'], 'file')
            disp(['Site ', num2str(s), ', Subject ', num2str(i) ': Computing  RISH features images for: ', data.ImageName]);
            
            %load dwi data
            [S, u, b, voxel, Sm, s0, ~, option] = loadDWI(data.ImageFullPath, ...
                data.ImageType);
            
            s0(isnan(s0))=0;
            s0(s0<=0)=1;
            
            
            %Remove possible noise
            S(S<0+eps) = 0;
            S(S>1-eps) = 1;
            S(isnan(S)) = 0;
            
            
            %load brain mask
            [mask, sd_mask, moption] = loadMask(data.MaskFullPath, data.MaskImageType, fid1);
            S = S.*single(repmat(mask,[1,1,1,size(S,4)]));
            
            % mask dwi signal and the average baseline
            s0 = s0.*single(mask);
            for j=1:size(S,4)
                S(:,:,:,j) = S(:,:,:,j).*s0;
            end
            
            newfile_suffix=[];
            if coption.denoise
                [S, Sigma] = MPdenoising(S, mask);
                newfile_suffix = [newfile_suffix, '_Denoised'];
                data.denoised = 1;
                
            end
            
            if max(b)~=coption.bValue && coption.bMap
                
                S = mapBValue(S, s0, b, coption.bValue);
                b = ones(size(S,4),1)*coption.bValue;
                
                newfile_suffix = [newfile_suffix, '_bMapped'];
                data.bValueMapped = 1;
                
            end
            
            newmask_suffix= ''; scale=[1 1 1];
            %rescale dwi data and brain mask
            if (abs(voxel(1)-coption.sp_high(1))>10^-3 ||abs(voxel(2)-coption.sp_high(2))>10^-3|| abs(voxel(3)-coption.sp_high(3))>10^-3) && coption.resample
                
                scale=coption.sp_high./voxel;
                [S, mask, voxel, Sm, s0, ~, sd_mask] = resampleData(S, mask, voxel, coption.sp_high,  Sm, s0, sd_mask, 7);
                
                newfile_suffix = [newfile_suffix, '_Resampled'];
                newmask_suffix= '_Resampled';
                data.resampled = 1;
                
            else
                
                s0(s0<=0) = 1;
                for j=1:size(S,4)
                    S(:,:,:,j) = S(:,:,:,j)./s0;
                    S(:,:,:,j) = S(:,:,:,j).*single(mask);
                end
                
            end
          
           
            Sout=repmat(s0,[1 1 1 size(S,4)]).*S;
            Sfinal = cat(4, s0, Sout);
            
            clear Sout;
           
            % if you applied any changes, save the new dwi and mask files
            if data.bValueMapped || data.resampled || data.denoised
                if strcmp(data.ImageType, 'nhdr') || strcmp(data.ImageType, 'nrrd')
                    disp('Saving new data...');
                    option.data=[];
                    dwirescaled = option;
                    dwirescaled.data = single(Sfinal);
                    dwirescaled.spacedirections = Sm(1:3,1:3);
                    dwirescaled.gradientdirections = [0 0 0; u(1:size(u,1)/2,:)];
                    dwirescaled.kinds = int32([1;1;1;4]);
                    if data.bValueMapped
                        dwirescaled.bvalue = coption.bValue;
                        dwirescaled.modality = 'DMRI';
                    end
                    nrrdSaveWithMetadata([data.Image_Harmonized_dir, data.ImageName, newfile_suffix, '.', data.ImageType],  dwirescaled);
                    
                    
                end
                
                
                if strcmp(data.MaskImageType, 'nhdr') || strcmp(data.MaskImageType, 'nrrd')
                    disp('Saving new mask...');
                    maskrescaled = moption;
                    maskrescaled.data = uint8(mask);
                    maskrescaled.spacedirections = sd_mask(1:3,1:3);
                    nrrdSaveWithMetadata([data.Image_Harmonized_dir, data.MaskName, newmask_suffix, '.', data.MaskImageType], maskrescaled);
                end
                
                
                if strcmp(data.ImageType, 'nii.gz') || strcmp(data.ImageType, 'nii')
                    disp('Saving new data...');
                    option.img=[];
                    dwirescaled = option;
                    dwirescaled.img = single(Sfinal);
                    dwirescaled.filetype=16;
                    dwirescaled.hdr.dime.datatype=16;
                    dwirescaled.hdr.dime.dim(5) = size(Sfinal,4);
                    dwirescaled.hdr.dime.dim(2:4)=[size(Sfinal,1) size(Sfinal,2) size(Sfinal,3)];
                    dwirescaled.hdr.dime.pixdim(2:4)=voxel(1:3);
                    
                    new_xyz = [option.hdr.hist.srow_x(1:3);option.hdr.hist.srow_y(1:3);option.hdr.hist.srow_z(1:3)]*diag(scale);
                    dwirescaled.hdr.hist.srow_x(1:3)=new_xyz(1,:);
                    dwirescaled.hdr.hist.srow_y(1:3)=new_xyz(2,:);
                    dwirescaled.hdr.hist.srow_z(1:3)=new_xyz(3,:);
                    
                    save_untouch_nii(dwirescaled, [data.Image_Harmonized_dir, data.ImageName, newfile_suffix, '.', data.ImageType]);
                    
                    fidBVec = fopen([data.Image_Harmonized_dir, data.ImageName, newfile_suffix, '.bvec'], 'w');
                    fprintf(fidBVec, '0 0 0\n');
                    for ui=1:size(u,1)/2
                        fprintf(fidBVec, '%f %f %f\n', u(ui,1), u(ui,2), u(ui,3));
                    end
                    fclose(fidBVec);
                    
                    fidBVal = fopen([data.Image_Harmonized_dir, data.ImageName, newfile_suffix, '.bval'], 'w');
                    fprintf(fidBVal, '0\n');
                    if coption.bMap
                        fprintf(fidBVal, '%f\n', repmat(coption.bValue,size(u,1)/2, 1) );
                    else
                        fprintf(fidBVal, '%f\n', repmat(max(b),size(u,1)/2, 1) );
                    end
                    fclose(fidBVal);
                    
                    
                end
                
                
                if strcmp(data.MaskImageType, 'nii.gz') || strcmp(data.MaskImageType, 'nii')
                    disp('Saving new mask...');
                    maskrescaled = moption;
                    maskrescaled.img = mask;
                    maskrescaled.hdr.dime.dim(2:4) = [size(S,1) size(S,2) size(S,3)];
                    maskrescaled.hdr.dime.pixdim(2:4) = voxel(1:3);%(diag(Sm(1:3,1:3)));
                    new_xyz = [moption.hdr.hist.srow_x(1:3);moption.hdr.hist.srow_y(1:3);moption.hdr.hist.srow_z(1:3)]*diag(scale);
                    maskrescaled.hdr.hist.srow_x(1:3) = new_xyz(1,:);
                    maskrescaled.hdr.hist.srow_y(1:3) = new_xyz(2,:);
                    maskrescaled.hdr.hist.srow_z(1:3) = new_xyz(3,:);
                    
                    save_untouch_nii(maskrescaled, [data.Image_Harmonized_dir, data.MaskName, newmask_suffix,'.', data.MaskImageType]);
                    
                    moption = maskrescaled;
                end
                
            end
            clear Sfinal;
            
            
            
            S = double(cat(4,S,S));
            
            
            
            %compute rish feature images and save them under Image_Harmonized_dir/rish/L0...L2...L4...L6...L8
            if(~isempty(S) && ~isempty(mask))
                disp('Computing RISH feature images...');
                % compute rish features
                IMG = rishFeatures(S, u,  mask, SHorder, coption.lambda);
                
                %save rish feature images
                if ~exist([data.Image_Harmonized_dir 'rish/'], 'dir')
                    mkdir([data.Image_Harmonized_dir 'rish/']);
                end
                
                
                % Convert rish feature images to nii.gz format (necessary
                % for ants registration)
                for l=1:SHorder/2+1
                    rishPath =  [data.Image_Harmonized_dir 'rish/', data.ImageName newfile_suffix '_L', num2str((l-1)*2)];
                    
                    
                    if strcmp(data.MaskImageType, 'nhdr') || strcmp(data.MaskImageType, 'nrrd') %nhdr/nrrd data type
                        
                        header = moption;
                        header.spacedirections = sd_mask(1:3,1:3);
                        header.data = single(IMG.rishImgs(:,:,:,l));
                        nrrdSaveWithMetadata([rishPath, '.nhdr'], header);
                        
                        % if your data is in .nhdr format, convert nhdr to nii since  ANTs based registration works only in nifti/analyze format
                        system(['ConvertBetweenFileFormats  ' rishPath  '.nhdr   ' rishPath '.nii.gz'  ]);
                        system(['rm  ' rishPath '.nhdr']);
                        system(['rm  ' rishPath '.raw.gz']);
                        
                    elseif strcmp(data.MaskImageType, 'nii.gz') || strcmp(data.MaskImageType, 'nii') %nifti type
                        moption.img = squeeze(IMG.rishImgs(:,:,:,l));
                        moption.filetype = 16;
                        moption.hdr.dime.datatype = 16;
                        moption.hdr.dime.dim(2:4) = [size(S,1) size(S,2) size(S,3)];
                        moption.hdr.dime.pixdim(2:4) = voxel(1:3);
                        new_xyz = [moption.hdr.hist.srow_x(1:3);moption.hdr.hist.srow_y(1:3);moption.hdr.hist.srow_z(1:3)]*diag(scale);
                        maskrescaled.hdr.hist.srow_x(1:3) = new_xyz(1,:);
                        maskrescaled.hdr.hist.srow_y(1:3) = new_xyz(2,:);
                        maskrescaled.hdr.hist.srow_z(1:3) = new_xyz(3,:);
                        save_untouch_nii(moption, [rishPath, '.nii.gz']);
                        
                    end
                    
                end
                nii = load_untouch_nii([rishPath '.nii.gz']); nii.img=s0.*double(mask);
                save_untouch_nii(nii,[data.Image_Harmonized_dir,  data.ImageName newfile_suffix '_b0.nii.gz']);
                
                
            end
            
            SITES{s}.InputImages{i} =  data; % update data (can be denoised/bvalue mapped/resampled)
            SITES{s}.InputImages{i}.DWISuffix  = newfile_suffix;
            SITES{s}.InputImages{i}.MaskSuffix = newmask_suffix;
            
            save([data.Image_Harmonized_dir,'newfile_suffix.mat'], 'newfile_suffix');
            save([data.Image_Harmonized_dir,'newmask_suffix.mat'], 'newmask_suffix');
            process = [data.denoised data.bValueMapped data.resampled];
            save([data.Image_Harmonized_dir,'processes.mat'], 'process');
            
            
            SITES{s}.InputImages{i}.denoised = process(1);
            SITES{s}.InputImages{i}.bValueMapped = process(2);
            SITES{s}.InputImages{i}.resampled = process(3);
        else
            
            load([data.Image_Harmonized_dir,'newfile_suffix.mat']);
            load([data.Image_Harmonized_dir,'newmask_suffix.mat']);
            load([data.Image_Harmonized_dir,'processes.mat']);
            SITES{s}.InputImages{i}.DWISuffix = newfile_suffix;
            SITES{s}.InputImages{i}.MaskSuffix = newmask_suffix;
            SITES{s}.InputImages{i}.denoised = process(1);
            SITES{s}.InputImages{i}.bValueMapped = process(2);
            SITES{s}.InputImages{i}.resampled = process(3);
            
        end
        
    end
end
fclose(fid1);
