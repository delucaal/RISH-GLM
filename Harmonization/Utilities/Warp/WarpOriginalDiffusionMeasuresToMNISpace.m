
%% This function harmonizes dMRI data for each site/scanner.
% Inputs:
% SITES: struct type of variable containing path information for each site and data
% templatesPath: the full path of the RISH feature templates
% Output:
% HarmonizedDWIPaths: contains the fullpath of each harmonized dMRI data

function  err=WarpOriginalDiffusionMeasuresToMNISpace(SITES, option)
err =0;
fid = fopen(['./log/log_', option.name,'WarpOriginalDiffusionMeasuresToMNISpace'],'w');

diffusionMeasures = {'FA','GFA', 'MD'};
for s=option.siteno:length(SITES)
    for i=1:SITES{s}.ImageNum
        data  = SITES{s}.InputImages{i}; % current data to be processed
        
        % data.DWISuffix: can be '_rescaled', '_bmapped', or '_denoised'
        fsldir = [data.ImageDirectory option.DTIdir  ];
        datapref = [data.ImageDirectory option.DTIdir data.ImageName  ];
        
        if ~exist([datapref  '_InMNI_FA.nii.gz'],'file') || option.force
            disp(['Registering Site: ', SITES{s}.name, ' Subject: ' num2str(i), ' ' data.ImageName, '...']);
            
            mov = [datapref, '_FA.nii.gz'];
            fix = './IITmean_FA.nii.gz';
            
            [status,out] = system(['antsRegistrationSyNQuick.sh -d 3 -f'...
                fix  ' -m '  mov  ' -o '  fsldir 'OrigToMNI_'  data.ImageName  ]);
            
            if status
                fprintf(['Error in antsRegistrationSyNQuick %s\n',  data.ImageFullPath]);
                fprintf('%s\n', out);
                err = 1; fclose(fid);
                return;
            end
            
            t_img1=[fsldir  'OrigToMNI_'  data.ImageName '1Warp.nii.gz'];
            t_img2=[fsldir  'OrigToMNI_'  data.ImageName  '0GenericAffine'];
            
            for mi=1:length(diffusionMeasures)
                
                mov_warped=[datapref, '_InMNI_', diffusionMeasures{mi}, '.nii.gz'];
                mov = [datapref, '_', diffusionMeasures{mi}, '.nii.gz'];
                
                [status,out] =system(['antsApplyTransforms -d 3 -i '  mov ' -o ' mov_warped ' -r  '  fix      '   -t '  t_img1  '    -t '  t_img2 '.mat' ]);
                
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

