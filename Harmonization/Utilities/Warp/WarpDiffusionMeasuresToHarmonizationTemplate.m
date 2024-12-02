
%% This function harmonizes dMRI data for each site/scanner.
% Inputs:
% SITES: struct type of variable containing path information for each site and data
% templatesPath: the full path of the RISH feature templates
% Output:
% HarmonizedDWIPaths: contains the fullpath of each harmonized dMRI data

function err= WarpDiffusionMeasuresToHarmonizationTemplate(SITES, templatesPath, option)
err =0;
fid = fopen(['./log/log_', option.name,'WarpDiffusionMeasuresToHarmonizationTemplate'],'w');
diffusionMeasures = {'FA','GFA', 'MD'};
for s=option.siteno:length(SITES)
    for i=1:SITES{s}.ImageNum
        data  = SITES{s}.InputImages{i}; % current data to be processed
        origName = data.ImageName;
        
        % data.DWISuffix: can be '_rescaled', '_bmapped', or '_denoised'
        if ~isempty(data.DWISuffix)
            data.ImageName = [data.ImageName, data.DWISuffix];
            data.MaskName = [data.MaskName, data.MaskSuffix];
            data.MaskFullPath =  [data.Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
            data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
            origName = data.ImageName;
        end
        
        if data.harmonized
            origName = data.ImageName; 
            data.ImageName = [option.harmonizedName, data.ImageName];
            data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
            data.MaskName = [option.harmonizedName, data.MaskName];
            data.MaskFullPath =  [data. Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
        end
        
        fsldir = [data.ImageDirectory option.DTIdir ];
        datapref = [data.ImageDirectory option.DTIdir data.ImageName  ];
        
        
        if ~exist([datapref  '_WarpedGFA.nii.gz'],'file') || option.force
            
            disp(['Registering Site: ', SITES{s}.name, ' Subject: ' num2str(i), ' ' data.ImageName, '...']);
            mov = [datapref, '_FA.nii.gz'];
            fix = [ templatesPath,    'Mean_' SITES{s}.name '_FA.nii.gz'];
            fixMask = [ templatesPath,     SITES{s}.name '_Mask.nii.gz'];
            if ~data.harmonized || ~exist([fsldir  'Pi_'  origName '1Warp.nii.gz'], 'file')
                [status,out] =system(['antsRegistrationSyNQuick.sh -d 3 -f'...
                    fix  ' -m '  mov ' -x ' fixMask ' -o '  fsldir 'Pi_'  origName  ]);
                if status
                    fprintf(['Error in antsRegistrationSyNQuick %s\n',  data.ImageFullPath]);
                    fprintf('%s\n', out);
                    err = 1; fclose(fid);
                    return;
                end
            end
        end
        
        if ~exist([datapref  '_WarpedFA.nii.gz'],'file') || ~exist([datapref  '_WarpedGFA.nii.gz'],'file') || ~exist([datapref  '_WarpedMD.nii.gz'],'file')
             %if ~option.training % if we are harmonizing the training data
                    fsldir = [data.ImageDirectory option.DTIdir ];
                    
                   
                    t_img1=[fsldir  'Pi_'  origName '1Warp.nii.gz'];
                    t_img2=[fsldir  'Pi_'  origName  '0GenericAffine.mat'];%MY017_NAA_024_L100GenericAffine;
%                 else
%                     [~, t_img1]=system(['ls ', templatesPath, origName '*1Warp.nii.gz']);
%                     [~, t_img2] = system(['ls ', templatesPath, origName '*0GenericAffine.mat']);
%                     t_img2=t_img2(1:end-1);  t_img1=t_img1(1:end-1);
%              end
%             
           
            fix = [ templatesPath,    'Mean_' SITES{s}.name '_FA.nii.gz'];
           
            for mi=1:length(diffusionMeasures)
                
                mov_warped=[datapref, '_Warped', diffusionMeasures{mi}, '.nii.gz'];
                mov = [datapref, '_', diffusionMeasures{mi}, '.nii.gz'];
                
                [status,out]=system(['antsApplyTransforms -d 3 -i '  mov ' -o ' mov_warped ' -r  '  fix      '   -t '  t_img1  '    -t '  t_img2  ]);
                
                if status
                    fprintf(['Error in antsApplyTransforms %s\n',  data.ImageFullPath]);
                    fprintf('%s\n', out);
                    err = 1; fclose(fid);
                    return;
                end
            end   
        end
    end
end
fclose(fid);
end

