
%% This function harmonizes dMRI data for each site/scanner.
% Inputs:
% SITES: struct type of variable containing path information for each site and data
% templatesPath: the full path of the RISH feature templates
% Output:
% err: contains the fullpath of each harmonized dMRI data

function  err=WarpHarmonizationTemplateToMNISpace(SITES, templatesPath, option)
err =0;
fid = fopen(['./log/log_', option.name,'WarpHarmonizationTemplateToMNISpace'],'w');

diffusionMeasures = {'FA','GFA', 'MD'};
for s=option.siteno:length(SITES)
    mov = [ templatesPath,    'Mean_' SITES{s}.name '_FA.nii.gz'];
    fix = './IITmean_FA.nii.gz';
   
    if ~exist([templatesPath  'TemplateToMNI_'  SITES{s}.name '1Warp.nii.gz'], 'file')
    system(['antsRegistrationSyNQuick.sh -d 3 -f'...
        fix  ' -m '  mov  ' -o '  templatesPath 'TemplateToMNI_'   SITES{s}.name  ]);
    end
    
    
    t_img1=[templatesPath  'TemplateToMNI_'  SITES{s}.name '1Warp.nii.gz'];
    t_img2=[templatesPath  'TemplateToMNI_'   SITES{s}.name  '0GenericAffine'];
    for i=1:SITES{s}.ImageNum
        
        data  = SITES{s}.InputImages{i}; % current data to be processed
        
       
        if ~isempty(data.DWISuffix)
            data.ImageName = [data.ImageName, data.DWISuffix];
            data.MaskName = [data.MaskName, data.MaskSuffix];
            data.MaskFullPath =  [data.Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
            data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
        end
        
        if data.harmonized
            data.ImageName = [option.harmonizedName, data.ImageName];
            data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
            data.MaskName = [option.harmonizedName, data.MaskName];
            data.MaskFullPath =  [data. Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
        end
      
        fsldir = [data.ImageDirectory option.DTIdir  ];
        datapref = [fsldir, data.ImageName  ];
        
        if ~exist([datapref  '_InMNI_FA.nii.gz'],'file') || option.force
            disp(['Registering Site: ', SITES{s}.name, ' Subject: ' num2str(i), ' ' data.ImageName, '...']);
            for mi=1:length(diffusionMeasures)
                
                mov_warped=[datapref, '_InMNI_', diffusionMeasures{mi}, '.nii.gz'];
                mov = [datapref, '_Warped', diffusionMeasures{mi}, '.nii.gz'];
                
                [status,out] = system(['antsApplyTransforms -d 3 -i '  mov ' -o ' mov_warped ' -r  '  fix      '   -t '  t_img1  '    -t '  t_img2 '.mat' ]);
                
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

