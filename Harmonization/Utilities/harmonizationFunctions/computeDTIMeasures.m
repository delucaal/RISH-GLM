% compute DTI measures using FSL and generalizedAnisotropy.m
function computeDTIMeasures(SITES, coption)

fid = fopen(['./log/log_', coption.name,'computeDTIMeasures'],'w');

for s=coption.siteno:length(SITES) % for all sites
    disp('');
    for i=1:SITES{s}.ImageNum % for all number of cases in each site
        
        data  = SITES{s}.InputImages{i}; % current data to be processed
        
        % data.DWISuffix: can be '_rescaled', '_bmapped', or '_denoised'
        if ~isempty(data.DWISuffix)
            data.ImageName = [data.ImageName, data.DWISuffix];
            data.MaskName = [data.MaskName, data.MaskSuffix];
            data.MaskFullPath =  [data. Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
            data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
        end
        
        if data.harmonized
            data.ImageName = [coption.harmonizedName, data.ImageName];
            data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
            data.MaskName = [coption.harmonizedName, data.MaskName];
            data.MaskFullPath =  [data.Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
        end
        
        fprintf(fid,'%s\n', data.ImageFullPath);
        %if data does not exist, exit
        
        % if data exists and data has not been processed before
        if ~exist([data.ImageDirectory coption.DTIdir data.ImageName  '_GFA.nii.gz'],'file') || ~exist([data.ImageDirectory coption.DTIdir data.ImageName  '_FA.nii.gz'],'file') || coption.force
            
            % create DTI measures
            disp(['Site ', num2str(s), ', Subject ', num2str(i) ': Creating DTI measures for: ', data.ImageName]);
            
            if ~exist([data.ImageDirectory coption.DTIdir], 'dir')
                mkdir([data.ImageDirectory coption.DTIdir]);
            end
            
            %% NRRD /NHDR data should be converted to nifti to use fsl
            if (strcmp(data.ImageType, 'nhdr') || strcmp(data.ImageType, 'nrrd'))
                
                % convert nrrd to nifti to use fsl dtifit
                [status,out] = system(['DWIConvert  --inputVolume ', data.ImageFullPath, ' -o ', data.ImageDirectory, 'temp.nii.gz', ...
                    ' --outputBVectors ', data.ImageDirectory, 'temp.bvec ',  ' --outputBValues ', data.ImageDirectory,...
                    'temp.bval --conversionMode NrrdToFSL --allowLossyConversion']);
                
                if status
                    fprintf(fid, '%s', out);  fclose(fid);
                    error(['Error in DWIConvert: During conversion of %s\n',  data.ImageFullPath]);
                end
                
                niitemp = load_untouch_nii([data.ImageDirectory, 'temp.nii.gz']);
                dwi = nrrdLoadWithMetadata(data.ImageFullPath);
                S =  dwi.data;
                u =  dwi.gradientdirections;
                switch size(u,1)
                    case size(S,4)
                        order = [1 2 3 4];
                    case size(S,1)
                        order = [2 3 4 1];
                    case size(S,2)
                        order = [1 4 3 2];
                    case size(S,3)
                        order = [1 2 4 3];
                end
                S = permute(S,order);
                
                niitemp.img =  S;
                save_untouch_nii(niitemp, [data.ImageDirectory, 'temp.nii.gz']);
                
                if (strcmp(data.MaskImageType, 'nhdr') || strcmp(data.MaskImageType, 'nrrd'))
                    iMask =  data.MaskFullPath; oMask = [ data.Image_Harmonized_dir, data.MaskName, '.nii.gz'];
                    [status,out] = system(['ConvertBetweenFileFormats  ',iMask, ' ',oMask]);
                    if status
                        fprintf(fid, ['Error in ConvertBetweenFileFormats: During conversion of %s\n',  out]);
                        fclose(fid);
                        error(['Error in ConvertBetweenFileFormats: %s\n',  data.ImageFullPath]);
                    end
                    
                    [status,out] = system(['module add fsl/6.0.1;dtifit -k ', data.ImageDirectory, 'temp.nii.gz -r ', data.ImageDirectory, 'temp.bvec -b ',...
                        data.ImageDirectory, 'temp.bval -m ', oMask, ' -o ' data.ImageDirectory coption.DTIdir, data.ImageName]);
                    if status
                        fprintf(fid, '%s', out); fclose(fid);
                        error(['Error in dtifit: %s\n',  data.ImageFullPath]);
                    end
                    
                    
                    
                    
                else
                    [status,out] = system(['dtifit -k ', data.ImageDirectory, 'temp.nii.gz -r ', data.ImageDirectory, 'temp.bvec -b ',...
                        data.ImageDirectory, 'temp.bval -m ', data.MaskFullPath, ' -o ' data.ImageDirectory coption.DTIdir, data.ImageName]);
                    if status
                        fprintf(fid, '%s', out); fclose(fid);
                        error(['Error in dtifit: %s\n',  data.ImageFullPath]);
                    end
                end
                
                 % remove mess
                system(['rm ', data.ImageDirectory, 'temp.nii.gz']);
                system(['rm ', data.ImageDirectory, 'temp.bval']);
                system(['rm ', data.ImageDirectory, 'temp.bvec']);
                %nifti type
            else
                if ~isempty(data.DWISuffix) || data.harmonized
                    [status,out] =system(['dtifit -k ', data.ImageFullPath, ' -m ', data.MaskFullPath,' -r ', data.Image_Harmonized_dir, data.ImageName '.bvec -b ',...
                        data.Image_Harmonized_dir, data.ImageName '.bval   -o ' data.ImageDirectory coption.DTIdir, data.ImageName]); 
                else
                     [status,out] =system(['dtifit -k ', data.ImageFullPath, ' -m ', data.MaskFullPath,' -r ', data.ImageDirectory, data.ImageName '.bvec -b ',...
                        data.ImageDirectory, data.ImageName '.bval   -o ' data.ImageDirectory coption.DTIdir, data.ImageName]);
                end
                  
               
                if status
                    fprintf(fid, '%s', out); fclose(fid);
                    error(['Error in dtifit: %s\n',  data.ImageFullPath]);
                end
            end
            
            % create GFA
            [S, ~, ~, ~, ~, ~, ~, ~] = loadDWI(data.ImageFullPath, ...
                data.ImageType);
            ga3D = generalizedAnisotropy(S);
            
            [mask, ~, ~] = loadMask(data.MaskFullPath, data.MaskImageType, 1);
            
            nii = load_untouch_nii([data.ImageDirectory coption.DTIdir, data.ImageName, '_FA.nii.gz']);
            nii.img = ga3D.*double(mask);
            save_untouch_nii(nii,[data.ImageDirectory, coption.DTIdir data.ImageName '_GFA.nii.gz']);
            
           
            
            
        end
        
        
    end
end
fclose(fid);











