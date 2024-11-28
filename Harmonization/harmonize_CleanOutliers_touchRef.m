%% This function harmonizes dMRI data for each site/scanner.
% Inputs:
% SITES: struct type of variable containing path information for each site and data
% templatesPath: the full path of the RISH feature templates
% Output:
% HarmonizedDWIPaths: contains the fullpath of each harmonized dMRI data

function [HarmonizedDWIPaths, SITES]  = harmonize_CleanOutliers_touchRef(SITES, templatesPath, coption)
fidH = fopen(['./log/log_', coption.name,'harmonization'],'w');

HarmonizedDWIPaths =[];

SHorder = coption.SHOrder; % Spherical Harmonics Order
for s=1:length(SITES) % harmonization is applied only to target  site. The other site gets rished
for do_harmonize=0:1 % when = 0 just produce the rished, when = 1 also harmonize
if(s == 1 && do_harmonize == 1)
    % not going to harmonize the ref
    continue
end
    
for i=1:SITES{s}.ImageNum
    data  = SITES{s}.InputImages{i}; % current data to be processed
    
    % data.DWISuffix: can be '_rescaled', '_bmapped', or '_denoised'
    if ~isempty(data.DWISuffix)
        data.ImageName = [data.ImageName, data.DWISuffix];
        data.MaskName = [data.MaskName, data.MaskSuffix];
        data.MaskFullPath =  [data.Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
        data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
    end
    if ~exist([data.Image_Harmonized_dir, coption.harmonizedName, data.ImageName,...
            '.', data.ImageType], 'file') || coption.force || strcmp(data.ImageName, 'GT_2669_dwi_xc_Ed_Resampled')
        
        
        mapped_Cs=[];
        
        disp(['Site:', SITES{s}.name, ' Subject ', num2str(i), ': ' data.ImageName]);
        imagepath = data.ImageFullPath;
        maskpath  = data.MaskFullPath;
        
        disp('Loading data...');
        [S, u, b, voxel, Sm, s0, ~, option]  = loadDWI(imagepath, ...
            data.ImageType);
        
        
        S = double(cat(4,S,S));
        S(S<0)=0; S(S>1)=1;
        s0(s0<=0)=1;
        
        disp('Loading mask...');
        [mask, ~, moption] = loadMask(maskpath, data.MaskImageType, 1);
        
        
        
        
        mov=[ templatesPath, 'Mean_', SITES{s}.name '_FA.nii.gz'];
        fixMask =maskpath;
        fix=[data.ImageDirectory coption.DTIdir   data.ImageName '_FA.nii.gz'];
        [status,out] =system(['antsRegistrationSyNQuick.sh -d 3 -f'...
            fix  ' -m '  mov ' -x ' fixMask ' -o '  data.Image_Harmonized_dir 'ToSubjectSpace_'  data.ImageName  ]);
        
        if status
            fprintf(fidH, '%s\n', out); fclose(fidH);
            error(['Error in antsRegistrationSyNQuick %s\n',  data.ImageFullPath]);
        end
        
        
        
        
        IMG = rishFeatures(S, u, mask, SHorder, coption.lambda); [sx,sy,sz,~]=size(IMG.rishImgs);
        
        rishPath = [data.Image_Harmonized_dir, 'rish/', data.ImageName];
        
        clear S;
        
        mapped_Cs = [];
        for L=0:2:SHorder
            
            disp(['Harmonizing L', num2str(L), '...']);
            mov=[templatesPath 'Scale_L' num2str(L) '.nii.gz'];
            
            
            
            fix=[ rishPath,    '_L'  num2str(L) '.nii.gz'];
            
            
            mov_warped=[data.Image_Harmonized_dir,     'rish/Scale_L' num2str(L) '_' data.ImageName  '.nii.gz'];
            t_img1=[data.Image_Harmonized_dir 'ToSubjectSpace_'  data.ImageName, '1Warp.nii.gz'];
            t_img2=[data.Image_Harmonized_dir 'ToSubjectSpace_'  data.ImageName, '0GenericAffine'];
            [status,out] =system(['antsApplyTransforms -d 3 -i '  mov ' -o ' mov_warped ' -r  '  fix      '   -t '  t_img1  '    -t '  t_img2 '.mat' ]);
            if status
                fprintf(fidH, '%s\n', out); fclose(fidH);
                error(['Error in antsApplyTransforms %s\n',  data.ImageFullPath]);
            end
                        
            
            
            Cs2_warped = load_untouch_nii(mov_warped);
            if(s == 1 || do_harmonize == 0)
               % ... don't scale!
               Cs2_warped.img = ones(size(Cs2_warped.img));
               disp('Not scaling');
            else
                disp('Scaling');
            end
%             disp('Warning!harmonize_cleanOutliers line 93');
%             Cs2_warped.img = sqrt(Cs2_warped.img);

            
            
            for SH_levels=[IMG.shs_same_level(L/2+1,1):IMG.shs_same_level(L/2+1,2)]
                mapped_Cs=[mapped_Cs ; (Cs2_warped.img(:)').*IMG.Cs(SH_levels,:)];
            end
        end
        
        S_hat= single(((IMG.Y*mapped_Cs))');
        S_hat_3D=reshape(S_hat,[sx sy sz length(u) ]);
        S_hat_3D(S_hat_3D<0) = 0;
        S_hat_3D(S_hat_3D>1) = 1;
    
        clear S_hat;
        if(s == 1 || do_harmonize == 0)
            IMG2 = IMG;
        else
            IMG2 = rishFeatures(S_hat_3D, u, mask, SHorder, 1e-5);%coption.lambda);
        end
        % Convert rish feature images to nii.gz format (necessary
        % for ants registration)
        disp('Saving harmonized RISH features...');
        for l=1:SHorder/2+1
            if(s == 1 || do_harmonize == 0)
                rishPathN =  [data.Image_Harmonized_dir 'rish/rished_', data.ImageName  '_L', num2str((l-1)*2)];
            else
                rishPathN =  [data.Image_Harmonized_dir 'rish/harmonized_', data.ImageName  '_L', num2str((l-1)*2)];
            end
            Cs2_warped.img = squeeze(IMG2.rishImgs(:,:,:,l));
            save_untouch_nii(Cs2_warped, [rishPathN, '.nii.gz']);
           
            % Sanity check
            t_img1=[data.Image_Harmonized_dir 'ToSubjectSpace_'  data.ImageName, '1InverseWarp.nii.gz'];
            t_img2=[data.Image_Harmonized_dir 'ToSubjectSpace_'  data.ImageName, '0GenericAffine'];
            [status,out] =system(['antsApplyTransforms -d 3 -i '  [rishPathN, '.nii.gz'] ' -o ' [rishPathN, '_Warped.nii.gz'] ' -r  '  mov      '   -t ['  t_img2  '.mat,1]    -t '  t_img1 '' ]);

        end
        
        S_hat_3D=repmat(s0,[1 1 1 size(S_hat_3D,4)]).*S_hat_3D;%
        Sfinal = cat(4, s0, S_hat_3D(:,:,:,1:size(S_hat_3D,4)/2));
        clear  S_hat_3D  IMG IMG2;
        disp('Saving harmonized dWI data...');
        % save harmonized dWI
        if(s == 1 || do_harmonize == 0)
            harmonizedDWI_img = [coption.rishedName data.ImageName];
        else
            harmonizedDWI_img = [coption.harmonizedName data.ImageName];
        end
        harmonizedDWI_dir = [data.Image_Harmonized_dir ];
        harmonizedDWI_type = data.ImageType;
        if strcmp(data.ImageType, 'nhdr') || strcmp(data.ImageType, 'nrrd')
            Output = option;
            Output.data = single(Sfinal);
            
            nrrdSaveWithMetadata([harmonizedDWI_dir,harmonizedDWI_img, '.', harmonizedDWI_type],  Output);
            
            
        elseif  strcmp(data.ImageType, 'nii.gz') || strcmp(data.ImageType, 'nii')
            Output = option;
            Output.hdr.dime.dim(5) = size(Sfinal,4);
            Output.img = single(Sfinal);
            
            save_untouch_nii(Output, [harmonizedDWI_dir,harmonizedDWI_img, '.', harmonizedDWI_type]);
            
            u = [0 0 0;u(1:end/2,:)];            
            u = u';
            save([harmonizedDWI_dir,harmonizedDWI_img, '.bvec'],'u','-ascii');
            b = [0 b'];
            save([harmonizedDWI_dir,harmonizedDWI_img, '.bval'],'b','-ascii');
            
            % ADL
%             fidBVec = fopen([harmonizedDWI_dir,harmonizedDWI_img, '.bvec'], 'w');
%             fprintf(fidBVec, '0 0 0\n');
%             for ui=1:size(u,1)/2
%                 fprintf(fidBVec, '%f %f %f\n', u(ui,1), u(ui,2), u(ui,3));
%             end
%             fclose(fidBVec);
            
%             fidBVal = fopen([harmonizedDWI_dir,harmonizedDWI_img, '.bval'], 'w');
%             fprintf(fidBVal, '0\n');
%             fprintf(fidBVal, '%f\n', repmat(coption.bValue,size(u,1)/2, 1) );
%             fclose(fidBVal);
            
            
        end
        clear Sfinal;
        SITES{s}.InputImages{i}.harmonized = 1;
        
        % outputs: Harmonized DWI paths and FA images
        
        HarmonizedDWIPaths{s, i} = [harmonizedDWI_dir   harmonizedDWI_img '.' harmonizedDWI_type];
        
        
    else
        SITES{s}.InputImages{i}.harmonized = 1;
    end
   
end
end
end
